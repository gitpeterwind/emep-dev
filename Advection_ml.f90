                    Module Advection_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DESCRIPTION
! 
! This module contains the routines for advection and diffusion.
! The sequence of advection in x y or z direction and vertical diffusion, 
! is controlled in the advecdiff routine.
!
! The horizontal advection is performed in the advx and advy routines using 
! "Bott's fourth order scheme". The routine preadvx and preadvy take care of 
! the transfer of information between processors before the advection step.
!
! The advvk routine performs the vertical advection. The advvdifvk does both
! vertical advection and diffusion. Bott's second order scheme with variable
! grid distance is used. The calculation of the coefficients used for this 
! scheme is done in the routine vgrid.
!
! Corrected version by pw, 1/11/01: xn_adv(k=top) is updated in 
! advvdifvk. We update fluxin and fluxout to account for the change in 
! concentrations which occur when the top layer is "reset".
! 
! Notes from Peter; 7/11/01
! About the division by the surface pressure (p*):
! In my opinion the best way is to divide by the advected p*, (corresponding 
! to the option where NSPEC_ADD_ADV is NOT equal to 0).  This should ensure 
! that, in the case of a uniform mixing ratio, we end up with a uniform mixing
! ratio, whatever the meteo.  The problem is that if the meteo is not 
! consistent (air is "created", or surface pressure does not corespond to the 
! quantity of air) the total weight of species may vary, creating problems for
! the mass budget. I see however no simple solution for this problem.
!
!
! Peter, January 2002: The advecdiff routine has been completely reorganised, 
! in order to allow for flexible timesteps. The timestep dt_advec can now be large. 
! The advecdiff routine will divide dt_advec in several advection steps if the 
! CFL condition is not met. For small dt_advec (600s in a 50x50 km2 grid)) 
! these changes should usually have no effect on the result.
! Some small changes have been made in the Bott's scheme in the advx, advy and 
! vgrid routines. The effect of these changes is estimated to be less than 1% on 
! average monthly values.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 use Par_ml            , only : MAXLIMAX,MAXLJMAX
 use ModelConstants_ml , only : KMAX_BND,KMAX_MID,NMET
!hf u2  use GenSpec_adv_ml, only: IXADV_O3,IXADV_H2 ,IXADV_NO2     ! DEBUG only
 implicit none
 private

 integer, public, parameter :: MIN_ADVGRIDS = 5 !minimum size of a subdomain

 integer, private, parameter :: NADVS      =  3  
 
 real, private, save, dimension(KMAX_BND)  ::  dhs1, dhs1i, dhs2i

!  for vertical advection (nonequidistant spacing)

 real, private, save, dimension(9,2:KMAX_MID,0:1)  ::  alfnew
 real, private, save, dimension(3)  ::  alfbegnew,alfendnew

 real, private,save, dimension(MAXLJMAX,KMAX_MID,NMET) :: uw,ue

 real, private,save, dimension(MAXLIMAX,KMAX_MID,NMET) :: vs,vn

 !ADVEC_TYPE (MADE=0, MACHO=1)
 integer, public, parameter :: ADVEC_TYPE = 1 ! Divides by p*

 public :: vgrid
 public :: advecdiff
 public :: adv_var
 public :: adv_int

 private :: advvk
 private :: advvdifvk
 private :: advx
 private :: advy
 private :: preadvx
 private :: preadvy

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	subroutine vgrid
!
!     inclusion of the variable grid spacing when interpolating the
!     polynominal is done by introducing new local coordinates, cor.
!
!
!     modified by pw january 2002: alfnew is modified such that 
!     a Courant number of one corresponds exactly to "empty" a cell. 
!     (small effects on results: less than 1%)
!


      use GridValues_ml, only : sigma_bnd,sigma_mid
	implicit none

	integer k,i,j
	real cor1, cor2, dcorl
      real hscor1(KMAX_BND+2),hscor2(KMAX_MID+2)
	real alfa1(NADVS), alfa2(NADVS), dei
      real corl1(KMAX_BND), corl2(KMAX_BND)

      real alf(9,2:KMAX_BND)

      do  k=1,KMAX_MID
        hscor1(k+1) = sigma_bnd(k)
        hscor2(k+1) = sigma_mid(k)
      enddo
      hscor1(KMAX_BND+1) = sigma_bnd(KMAX_BND)

      hscor1(1) = - sigma_bnd(2)
      hscor1(KMAX_BND+2) = 2.*sigma_bnd(KMAX_BND) - sigma_bnd(KMAX_BND-1)
      hscor2(1) = 2.*sigma_mid(1) - sigma_mid(2)
      hscor2(KMAX_MID+2) = 2.*sigma_mid(KMAX_MID) - sigma_mid(KMAX_MID-1)     

      do  k=1,KMAX_BND
	  dhs1(k) = hscor1(k+1) - hscor1(k)
	  dhs1i(k) = 1./dhs1(k)
	  dhs2i(k) = 1./(hscor2(k+1) - hscor2(k))
      enddo

      do k=2,KMAX_BND

	  corl1(k) = (hscor1(k) - hscor2(k))*dhs1i(k)
	  corl2(k) = (hscor1(k+1) - hscor2(k))*dhs1i(k)
	  dcorl = corl2(k) - corl1(k)
	  alfa1(NADVS-1) = (corl2(k)**2 - corl1(k)**2)/(2.*dcorl)
	  alfa2(NADVS-1) = (corl2(k)**3 - corl1(k)**3)/(3.*dcorl)

	  cor1 = (hscor1(k-1) - hscor2(k))*dhs1i(k)
!	  cor2 = (hscor1(k) - hscor2(k))*dhs1i(k)
	  dcorl = corl1(k) - cor1
	  alfa1(NADVS-2) = (corl1(k)**2 - cor1**2)/(2.*dcorl)
	  alfa2(NADVS-2) = (corl1(k)**3 - cor1**3)/(3.*dcorl)
     
!	  cor1 = (hscor1(k+1) - hscor2(k))*dhs1i(k)
	  cor2 = (hscor1(k+2) - hscor2(k))*dhs1i(k)
	  dcorl = cor2 - corl2(k)
	  alfa1(NADVS) = (cor2**2 - corl2(k)**2)/(2.*dcorl)
	  alfa2(NADVS) = (cor2**3 - corl2(k)**3)/(3.*dcorl)

	  dei = alfa1(NADVS-1)*alfa2(NADVS) - alfa1(NADVS)*alfa2(NADVS-1) &
		- alfa1(NADVS-2)*(alfa2(NADVS) - alfa2(NADVS-1))	  &
		+ alfa2(NADVS-2)*(alfa1(NADVS) - alfa1(NADVS-1))
	  dei = 1./dei

	  alf(1,k) = dei	&
		*(alfa1(NADVS-1)*alfa2(NADVS)-alfa1(NADVS)*alfa2(NADVS-1))
	  alf(4,k) = dei	&
		*(alfa2(NADVS-2)*alfa1(NADVS)- alfa2(NADVS)*alfa1(NADVS-2))
	  alf(7,k) = dei	&
		*(alfa2(NADVS-1)*alfa1(NADVS-2)-alfa2(NADVS-2)*alfa1(NADVS-1))
	  alf(2,k) = dei	&
		*(alfa2(NADVS-1) - alfa2(NADVS))/2.
	  alf(5,k) = dei	&
		*(alfa2(NADVS) - alfa2(NADVS-2))/2.
	  alf(8,k) = dei	&
		*(alfa2(NADVS-2) - alfa2(NADVS-1))/2.
	  alf(3,k) = dei	&
		*(alfa1(NADVS) - alfa1(NADVS-1))/3.
	  alf(6,k) = dei	&
		*(alfa1(NADVS-2) - alfa1(NADVS))/3.
	  alf(9,k) = dei	&
		*(alfa1(NADVS-1) - alfa1(NADVS-2))/3.
	enddo

      do k=2,KMAX_MID
	  alfnew(1,k,0) = alf(1,k) + 2.*alf(2,k)*corl2(k)		&
		+ 3.*alf(3,k)*corl2(k)*corl2(k)
	  alfnew(1,k,1) = -(alf(1,k+1) + 2.*alf(2,k+1)*corl1(k+1)	&
		+ 3.*alf(3,k+1)*corl1(k+1)*corl1(k+1))
	  alfnew(4,k,0) = alf(4,k) + 2.*alf(5,k)*corl2(k)		&
		+ 3.*alf(6,k)*corl2(k)*corl2(k)
	  alfnew(4,k,1) = -(alf(4,k+1) + 2.*alf(5,k+1)*corl1(k+1)	&
		+ 3.*alf(6,k+1)*corl1(k+1)*corl1(k+1))
	  alfnew(7,k,0) = alf(7,k) + 2.*alf(8,k)*corl2(k)		&
		+ 3.*alf(9,k)*corl2(k)*corl2(k)
	  alfnew(7,k,1) = -(alf(7,k+1) + 2.*alf(8,k+1)*corl1(k+1)	&
		+ 3.*alf(9,k+1)*corl1(k+1)*corl1(k+1))
!pw	  alfnew(2,k,0) = -(alf(2,k) + 3.*alf(3,k)*corl2(k))		&
!		*dhs2i(k)
!	  alfnew(2,k,1) = (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))	&
!		*dhs2i(k)
	  alfnew(2,k,0) = -(alf(2,k) + 3.*alf(3,k)*corl2(k))		&
		*dhs1i(k)
	  alfnew(2,k,1) = (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))	&
		*dhs1i(k+1)
!pw	  alfnew(5,k,0) = -(alf(5,k) + 3.*alf(6,k)*corl2(k))		&
!		*dhs2i(k)
!	  alfnew(5,k,1) = (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))	&
!		*dhs2i(k)
	  alfnew(5,k,0) = -(alf(5,k) + 3.*alf(6,k)*corl2(k))		&
		*dhs1i(k)
	  alfnew(5,k,1) = (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))	&
		*dhs1i(k+1)
!pw	  alfnew(8,k,0) = -(alf(8,k) + 3.*alf(9,k)*corl2(k))		&
!		*dhs2i(k)
!	  alfnew(8,k,1) = (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))	&
!		*dhs2i(k)
	  alfnew(8,k,0) = -(alf(8,k) + 3.*alf(9,k)*corl2(k))		&
		*dhs1i(k)
	  alfnew(8,k,1) = (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))	&
		*dhs1i(k+1)
!pw	  alfnew(3,k,0) = alf(3,k)*dhs2i(k)*dhs2i(k)
!	  alfnew(3,k,1) = -alf(3,k+1)*dhs2i(k)*dhs2i(k)
!	  alfnew(6,k,0) = alf(6,k)*dhs2i(k)*dhs2i(k)
!	  alfnew(6,k,1) = -alf(6,k+1)*dhs2i(k)*dhs2i(k)
!	  alfnew(9,k,0) = alf(9,k)*dhs2i(k)*dhs2i(k)
!	  alfnew(9,k,1) = -alf(9,k+1)*dhs2i(k)*dhs2i(k)
	  alfnew(3,k,0) = alf(3,k)*dhs1i(k)*dhs1i(k)
	  alfnew(3,k,1) = -alf(3,k+1)*dhs1i(k+1)*dhs1i(k+1)
	  alfnew(6,k,0) = alf(6,k)*dhs1i(k)*dhs1i(k)
	  alfnew(6,k,1) = -alf(6,k+1)*dhs1i(k+1)*dhs1i(k+1)
	  alfnew(9,k,0) = alf(9,k)*dhs1i(k)*dhs1i(k)
	  alfnew(9,k,1) = -alf(9,k+1)*dhs1i(k+1)*dhs1i(k+1)
	enddo

	alfbegnew(1) = alfnew(1,2,0)+alfnew(4,2,0)
	alfbegnew(2) = alfnew(2,2,0)+alfnew(5,2,0)
	alfbegnew(3) = alfnew(3,2,0)+alfnew(6,2,0)
      alfendnew(1) = alfnew(4,KMAX_MID,1)+alfnew(7,KMAX_MID,1)
      alfendnew(2) = alfnew(5,KMAX_MID,1)+alfnew(8,KMAX_MID,1)
      alfendnew(3) = alfnew(6,KMAX_MID,1)+alfnew(9,KMAX_MID,1)

	end subroutine vgrid

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	subroutine advecdiff
!___________________________________________________________________________________
!Flexible timestep. Peter Wind january-2002
!
! dt_advec : time interval between two advections calls 
! (controls time splitting between advection and chemistry )
!
! dt_xys : time intervall between vertical and horizontal advection steps
!(controls time splitting between vertical and horizontal advection)
! There is one sequence (z),(x,y,y,x),(z) during each dt_xys 
!
! dt_xy : time intervall for horizontal advection iterations
!(controls time splitting between x and y advection)
! There is one sequence x,y,y,x during each dt_xy
!
! dt_s : time intervall for vertical advection iterations
!
! dt_advec >= dt_xys >= max(dt_xy, dt_s)
!
	use Par_ml   , only : li0,li1,lj0,lj1	&
			,limax,ljmax		&
			,me,NPROC
	use GenSpec_adv_ml , only : NSPEC_ADV
        use ModelConstants_ml, only : nstep, dt_advec, 	&
			PT,KCHEMTOP
        use GridValues_ml, only : GRIDWIDTH_M,xm2,xmd,xm2ji,xmdji,carea 
	use Met_ml ,only : ps,sdot,skh,u,v
	use Chemfields_ml, only : xn_adv
	use My_Timing_ml,  only : Code_timer, Add_2timing, tim_before, tim_after
        use MassBudget_ml , only : fluxin,fluxout

	implicit none

!	local

	integer i,j,k,n,l,info
	real dtsave,dth
	real xntop(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
	real xnw(3*NSPEC_ADV,MAXLJMAX),xne(3*NSPEC_ADV,MAXLJMAX)
	real xnn(3*NSPEC_ADV,MAXLIMAX), xns(3*NSPEC_ADV,MAXLIMAX)
      real ps3d(MAXLIMAX,MAXLJMAX,KMAX_MID)
	real psw(3,MAXLJMAX),pse(3,MAXLJMAX)
	real psn(3,MAXLIMAX), pss(3,MAXLIMAX)
      real ds3(2:KMAX_MID),ds4(2:KMAX_MID)
	real ulmin,ulmax,vlmin,vlmax,ucmax,vcmax
	integer inadvst
	logical lvertsplit
	real dhskmax,sdotmax,sdotmin
	real sdotmaxk,sdotmink
	real sdotmaxadv,sum1

        real xcmax(KMAX_MID),ycmax(KMAX_MID),scmax,sdcmax
        real dt_xysmax,dt_xymax(KMAX_MID),dt_smax
        real dt_xys,dt_xy(KMAX_MID),dt_s
        integer niterxys,niterxy(KMAX_MID),niters,nxy
        integer iterxys,iterxy,iters

        call Code_timer(tim_before)

	if(KCHEMTOP==2)then
	  xntop(:,:,:)=xn_adv(:,:,:,1)
	endif

!   convert from mixing ratio to concentration before advection

      do k = 1,KMAX_MID
	  do j = 1,ljmax
	    do i = 1,limax

	      xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*(ps(i,j,1)-PT)

	      ps3d(i,j,k) = ps(i,j,1) - PT

	    end do
	  end do
	end do

        call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

!     time-splitting is used for the physical and chemical operators.
!     second-order accuracy in time is obtained by alternating the order 
!     of the advx and advy operators from one time-step to another.

!
! Determine timestep for horizontal advection.
! 
! Courant criterion, which takes into account the mapping factor xm2:
! left face:    xm2(i)    u(i)  dt/dx < 1  when u(i) > 0
! right face:   xm2(i) |u(i-1)| dt/dx < 1  when u(i-1) < 0
!
! In the case where the flux is streaming out of the cell i from both faces, 
! then the total should be < 1:
! xm2(i) |u(i-1)| dt/dx + xm2(i) u(i) dt/dx < 1    for u(i-1)<0 and u(i)>0
!
! The three conditions can be written as:
!
! max(xm2(i)*u(i)*dt/dx , 0.0) - min(xm2(i)*u(i-1)*dt/dx , 0.0) < 1
!
! or equivalently:
! dt < dx / ( max(xm2(i)*u(i) , 0.0) - min(xm2(i)*u(i-1) , 0.0) )
!
! In the case of variable cell size, dx is defined as dx(i) in these formula.
!
! The value 1.e-30 is to ensure that we don't divide by 0 when all the velocities are 0.

	dth = dt_advec/GRIDWIDTH_M

        do k=1,KMAX_MID

           xcmax(k) = maxval(amax1(u(1:limax,1:ljmax,k,1)*xm2(1:limax,1:ljmax), 1.e-30) &
                          -amin1(u(0:limax-1,1:ljmax,k,1)*xm2(1:limax,1:ljmax), 0.0) )
           ycmax(k) = maxval(amax1(v(1:limax,1:ljmax,k,1)*xm2(1:limax,1:ljmax), 0.0) &
                          -amin1(v(1:limax,0:ljmax-1,k,1)*xm2(1:limax,1:ljmax), 0.0) )
           xcmax(k) = amax1(xcmax(k),ycmax(k))
        enddo
	  call gc_rmax(KMAX_MID, NPROC, info, xcmax)
          

        do k=1,KMAX_MID
          dt_xymax(k)=GRIDWIDTH_M/xcmax(k)
        enddo
44      format('k =',I4,6F12.4)



!Courant number in vertical sigma coordinates:  sigmadot*dt/deltasigma
!
!Note that dhs1(k+1) denotes thickness of layer k 
!     and sdot(k+1) denotes sdot at the boundary between layer k and k+1
!
!flux through wall k+1:  sdot(k+1) *dt/dhs1(k+1)<1   for sdot(k+1)>0 
!                       |sdot(k+1)|*dt/dhs1(k+2)<1   for sdot(k+1)<0 
!
!layer k: sdot(k+1)*dt/dhs1(k+1) + |sdot(k)|*dt/dhs1(k+1) <1 for sdot(k+1)>0 and sdot(k)<0 
! 
!total out of layer k: amax1(sdot(1:limax,1:ljmax,k+1,1),0.0)-amin1(sdot(1:limax,1:ljmax,k,1),0.0) 
!
        scmax = 1.e-30
        do k = 1,KMAX_MID
           sdcmax=maxval( amax1(sdot(1:limax,1:ljmax,k+1,1),0.0) &
                -amin1(sdot(1:limax,1:ljmax,k,1),0.0)   )
           scmax = amax1(sdcmax/dhs1(k+1),scmax)
        enddo

        call gc_rmax(1, NPROC, info, scmax)
        dt_smax = 1./scmax
           
        dt_xysmax = amax1(dt_smax, maxval(dt_xymax(1:KMAX_MID)))
        niterxys = int(dt_advec/dt_xysmax)+1
        dt_xys = dt_advec/real(niterxys)
        niters = int(dt_xys/dt_smax)+1
        dt_s = dt_xys/real(niters)
!        if(me.eq.0)then
!           write(*,45)dt_xysmax,dt_xys,niterxys
!           write(*,45)dt_smax,dt_s,niters
!        endif
        nxy=0
        do k=1,KMAX_MID
           niterxy(k) = int(dt_xys/dt_xymax(k))+1
           dt_xy(k) = dt_xys/real(niterxy(k))
           nxy=nxy+niterxy(k)
!        if(me.eq.0)then
!        write(*,46)k,dt_xymax(k),dt_xy(k),niterxy(k)
!        endif
        enddo
        if(me.eq.0)then
        write(*,47)niterxys-1,nxy-KMAX_MID,niters-1
        endif
!        if(me.eq.0)then
!        write(*,47)nxy,nxy/20.
!        endif
45      format(2F12.4,I6)
46      format('k = ',I6,2F12.4,I6,2F12.4,I6)
47      format('extra iterations (xyz,xy,z): ',3I6)

        call Add_2timing(20,tim_after,tim_before,	&
			"advecdiff:initialisations")

! Start xys advection loop:
        iterxys = 0
        do while (iterxys < niterxys) 
           if(mod(nstep,2) /= 0 .or. iterxys /= 0)then !start a xys sequence
              
              iterxys = iterxys + 1
              do k = 1,KMAX_MID
                 dth = dt_xy(k)/GRIDWIDTH_M
              do iterxy=1,niterxy(k)

                 ! send/receive in x-direction

                 call preadvx(1100+k				&
                      ,xn_adv(1,1,1,k),ps3d(1,1,k),u(0,1,k,1)	&
                      ,xnw,xne				&
                      ,psw,pse)

                 !	x-direction
                 do j = lj0,lj1
                    call advx(					&
                         u(0,j,k,1),uw(j,k,1),ue(j,k,1)		&
                         ,xn_adv(1,1,j,k),xnw(1,j),xne(1,j)	&
                         ,ps3d(1,j,k),psw(1,j),pse(1,j)		&
                         ,xm2(0,j),xmd(0,j)			&
                         ,dth,carea(k))

                 enddo

                 call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx")

                 !	send/receive in y-direction

                 call preadvy(1300+k				&
                      ,xn_adv(1,1,1,k),ps3d(1,1,k),v(1,0,k,1)	&
                      ,xns, xnn				&
                      ,pss, psn)

                 !	y-direction

                 call Add_2timing(22,tim_after,tim_before,		&
                      "advecdiff:preadvy")

                 do i = li0,li1
                    call advy(					&
                         v(i,0,k,1),vs(i,k,1),vn(i,k,1)		&
                         ,xn_adv(1,i,1,k),xns(1,i),xnn(1,i)	&
                         ,ps3d(i,1,k),pss(1,i),psn(1,i)		&
                         ,xm2ji(0,i),xmdji(0,i)			&
                         ,dth,carea(k))

                 enddo

                 call Add_2timing(23,tim_after,tim_before,"advecdiff:advy")

              enddo !iterxy horizontal (xy) advection
              enddo !k horizontal (xy) advection

              do iters=1,niters
                 if(iterxys < niterxys .or. iters < niters )then ! The last vertical advection call is done separately

                    !perform only vertical advection
                    do j = lj0,lj1
                       do i = li0,li1

                          call advvk(xn_adv(1,i,j,1),ps3d(i,j,1)	&
                               ,sdot(i,j,1,1),dt_s)

                       enddo
                    enddo

                 endif ! 

              enddo ! vertical (s) advection
                      call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

           else  !start a yxs sequence

              iterxys = iterxys + 1

              do k = 1,KMAX_MID
                 dth = dt_xy(k)/GRIDWIDTH_M
              do iterxy=1,niterxy(k)

                 !	send/receive in y-direction

                 call preadvy(1300+k				&
                      ,xn_adv(1,1,1,k),ps3d(1,1,k),v(1,0,k,1)	&
                      ,xns, xnn				&
                      ,pss, psn)

                 !	y-direction

                 call Add_2timing(22,tim_after,tim_before,		&
                      "advecdiff:preadvy")

                 do i = li0,li1
                    call advy(					&
                         v(i,0,k,1),vs(i,k,1),vn(i,k,1)		&
                         ,xn_adv(1,i,1,k),xns(1,i),xnn(1,i)	&
                         ,ps3d(i,1,k),pss(1,i),psn(1,i)		&
                         ,xm2ji(0,i),xmdji(0,i)			&
                         ,dth,carea(k))

                 enddo

                 call Add_2timing(23,tim_after,tim_before,"advecdiff:advy")

                 ! send/receive in x-direction

                 call preadvx(1100+k				&
                      ,xn_adv(1,1,1,k),ps3d(1,1,k),u(0,1,k,1)	&
                      ,xnw,xne				&
                      ,psw,pse)

                 !	x-direction
                 do j = lj0,lj1
                    call advx(					&
                         u(0,j,k,1),uw(j,k,1),ue(j,k,1)		&
                         ,xn_adv(1,1,j,k),xnw(1,j),xne(1,j)	&
                         ,ps3d(1,j,k),psw(1,j),pse(1,j)		&
                         ,xm2(0,j),xmd(0,j)			&
                         ,dth,carea(k))

                 enddo

                 call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx")

              enddo !iterxy horizontal (xy) advection
              enddo !k horizontal (xy) advection

              do iters=1,niters
                 if(iterxys < niterxys .or. iters < niters )then ! The last vertical advection call is done separately

                    !perform only vertical advection
                    do j = lj0,lj1
                       do i = li0,li1

                          call advvk(xn_adv(1,i,j,1),ps3d(i,j,1)	&
                               ,sdot(i,j,1,1),dt_s)

                       enddo
                    enddo

                 endif ! 

              enddo ! vertical (s) advection
                      call Add_2timing(24,tim_after,tim_before,"advecdiff:advvkdiff")

           endif ! yxs sequence
        enddo

        !last vertical advection call. Perform also division by p* and diffusion.

        do k = 2,KMAX_MID
           ds3(k) = dt_advec*dhs1i(k)*dhs2i(k)
           ds4(k) = dt_advec*dhs1i(k+1)*dhs2i(k)
        enddo

        do j = lj0,lj1
           do i = li0,li1
              call advvdifvk(xn_adv(1,i,j,1),ps3d(i,j,1)	&
                   ,sdot(i,j,1,1),skh(i,j,1,1)		&
                   ,ds3,ds4				&
                   ,ps(i,j,1)-PT,dt_s)

           enddo
        enddo

	if(lj0.ne.1)then 
        do k=KCHEMTOP,KMAX_MID
	    do i = 1,limax
		xn_adv(:,i,1,k) = xn_adv(:,i,1,k)/(ps(i,1,1)-PT)
	    enddo
	  enddo
	endif
	if(li0.ne.1)then
        do k=KCHEMTOP,KMAX_MID
            do j=lj0,lj1    
		xn_adv(:,1,j,k) = xn_adv(:,1,j,k)/(ps(1,j,1)-PT)
            enddo
	  enddo
	endif         
	if(li1.ne.limax)then
        do k=KCHEMTOP,KMAX_MID
            do j=lj0,lj1    
		xn_adv(:,limax,j,k) = xn_adv(:,limax,j,k)	&
			/(ps(limax,j,1)-PT)
            enddo
	  enddo
	endif         
	if(lj1.ne.ljmax)then 
        do k=KCHEMTOP,KMAX_MID
            do i = 1,limax
		xn_adv(:,i,ljmax,k) = xn_adv(:,i,ljmax,k)	&
			/(ps(i,ljmax,1)-PT)
            enddo
	  enddo
	endif 
             
	if(KCHEMTOP==2)then

!pw since the xn_adv are changed it corresponds to a flux in or 
!   out of the system:  
 
           do i  = li0,li1
              do j = lj0,lj1
           where(xn_adv(:,i,j,1) .gt. xntop(:,i,j))
              fluxout(:) = fluxout(:) + &
                   (xn_adv(:,i,j,1) - xntop(:,i,j)) &
                   *(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
           else where
              fluxin(:) = fluxin(:) + &
                   (xntop(:,i,j) - xn_adv(:,i,j,1)) &
                   *(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
           end where


              enddo
           enddo
	  xn_adv(:,:,:,1) = xntop(:,:,:)

	endif

        call Add_2timing(24,tim_after,tim_before,"advecdiff:advvkdiff")


	end subroutine advecdiff

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine advvk(xn_adv,ps3d,sdot,dt_s)

!     executes advection with a. bott's integreated flux-form
!     using 2'nd order polynomial in the vertical.

	use ModelConstants_ml   , only : EPSIL, dt_advec
	use GenSpec_adv_ml , only : NSPEC_ADV
	implicit none

!	input
      real,intent(in)::  sdot(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1),dt_s

!	input+output
      real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)
      real ,intent(inout):: ps3d(0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)

!	local
      real fluxk(NSPEC_ADV,KMAX_MID),fluxps(KMAX_MID),fc(KMAX_MID)

!	local
	integer  k, n1k,k1
	integer klimlow,klimhig
	real zzfl1,zzfl2,zzfl3,totk(NSPEC_ADV),totps
	real fc1,fc2,fc3

      do k = 1,KMAX_MID-1
	  fc(k) = sdot(k*MAXLIMAX*MAXLJMAX)*dt_s
	enddo

      fc(KMAX_MID) = -1.
!-------------- calculate the advection ----------------------------

	klimlow = 1
	if(fc(1).ge.0.)klimlow=2
      klimhig = KMAX_MID-1
      if(fc(KMAX_MID-1).lt.0.)klimhig = KMAX_MID-2

	fluxk(:,1) = 0.
	fluxps(1) = 0.

	if(fc(1).ge.0.)then

	  fc1 = fc(1)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
	  zzfl2 = alfbegnew(1)*fc1		&
		 + alfbegnew(2)*fc2		&
		 + alfbegnew(3)*fc3
	  zzfl3 = alfnew(7,2,0)*fc1		&
		 + alfnew(8,2,0)*fc2		&
		 + alfnew(9,2,0)*fc3

	  fluxk(:,2) = amax1(0.,xn_adv(:,0)*zzfl2	&
		+xn_adv(:,MAXLIMAX*MAXLJMAX)*zzfl3)
	  fluxps(2) = amax1(0.,ps3d(0)*zzfl2		&
		+ps3d(MAXLIMAX*MAXLJMAX)*zzfl3)

	endif
	do k = klimlow,klimhig

	  fc1 = fc(k)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
	  n1k = 0
	  if(fc1.lt.0)n1k=1
	  zzfl1 = alfnew(1,k,n1k)*fc1		&
		 + alfnew(2,k,n1k)*fc2		&
		 + alfnew(3,k,n1k)*fc3
	  zzfl2 = alfnew(4,k,n1k)*fc1		&
		 + alfnew(5,k,n1k)*fc2		&
		 + alfnew(6,k,n1k)*fc3
	  zzfl3 = alfnew(7,k,n1k)*fc1		&
		 + alfnew(8,k,n1k)*fc2		&
		 + alfnew(9,k,n1k)*fc3

	  k1 = k-1+n1k

	  fluxk(:,k+1) = amax1(0.,xn_adv(:,(k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1 &
		+xn_adv(:,k1*MAXLIMAX*MAXLJMAX)*zzfl2			 &
		+xn_adv(:,(k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)
	  fluxps(k+1) = amax1(0.,ps3d((k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1	&
		+ps3d(k1*MAXLIMAX*MAXLJMAX)*zzfl2			&
		+ps3d((k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)

	enddo
      if(fc(KMAX_MID-1).lt.0.)then

        fc1 = fc(KMAX_MID-1)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
        zzfl1 = alfnew(1,KMAX_MID,1)*fc1            &
             + alfnew(2,KMAX_MID,1)*fc2      &
             + alfnew(3,KMAX_MID,1)*fc3
	  zzfl2 = alfendnew(1)*fc1		&
		 + alfendnew(2)*fc2		&
		 + alfendnew(3)*fc3

        fluxk(:,KMAX_MID) =                                     &
            amax1(0.,xn_adv(:,(KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1      &
            +xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2) 
        fluxps(KMAX_MID) =                                     &
            amax1(0.,ps3d((KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1      &
            +ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2) 

	endif

	k=1

      do while(k.lt.KMAX_MID)

	  if(fc(k).lt.0.) then

	    if(fc(k+1).ge.0.) then
	        totk(:) = amin1(xn_adv(:,k*MAXLIMAX*MAXLJMAX)	&
			*dhs1(k+2)/					&
	        	(fluxk(:,k+1) + fluxk(:,k+2)+ EPSIL),1.)
	        fluxk(:,k+1) = -fluxk(:,k+1)*totk(:)
	        fluxk(:,k+2) = fluxk(:,k+2)*totk(:)
	        xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) = 			&
			amax1(0., xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)	&
	                -(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
	        xn_adv(:,k*MAXLIMAX*MAXLJMAX) = 			&
			amax1(0., xn_adv(:,k*MAXLIMAX*MAXLJMAX)	&
	                 -(fluxk(:,k+2) - fluxk(:,k+1))*dhs1i(k+2))

	        totps = amin1(ps3d(k*MAXLIMAX*MAXLJMAX)		&
			*dhs1(k+2)/					&
			(fluxps(k+1) + fluxps(k+2)+ EPSIL),1.)
	        fluxps(k+1) = -fluxps(k+1)*totps
	        fluxps(k+2) = fluxps(k+2)*totps
	        ps3d((k-1)*MAXLIMAX*MAXLJMAX) = 			&
			amax1(0., ps3d((k-1)*MAXLIMAX*MAXLJMAX)		&
	                -(fluxps(k+1) - fluxps(k))*dhs1i(k+1))
	        ps3d(k*MAXLIMAX*MAXLJMAX) = 			&
			amax1(0., ps3d(k*MAXLIMAX*MAXLJMAX)		&
			-(fluxps(k+2) - fluxps(k+1))*dhs1i(k+2))
	      k = k+2
	    else
	      fluxk(:,k+1) = 						&
		-amin1(xn_adv(:,k*MAXLIMAX*MAXLJMAX)*dhs1(k+2),	&
		       	fluxk(:,k+1))
	      xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) = 			&
		amax1(0., xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)		&
			-(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
	      fluxps(k+1) = 						&
		-amin1(ps3d(k*MAXLIMAX*MAXLJMAX)*dhs1(k+2),		&
			fluxps(k+1))
	      ps3d((k-1)*MAXLIMAX*MAXLJMAX) = 				&
		amax1(0., ps3d((k-1)*MAXLIMAX*MAXLJMAX)			&
			-(fluxps(k+1) - fluxps(k))*dhs1i(k+1))

	      k = k+1
	    endif
	  else

	    fluxk(:,k+1) = amin1(xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1), &
			fluxk(:,k+1))
	    xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) = 			&
		amax1(0., xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)		&
	                 -(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
	    fluxps(k+1) = amin1(ps3d((k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1),	&
		fluxps(k+1))
	    ps3d((k-1)*MAXLIMAX*MAXLJMAX) = 				&
		amax1(0., ps3d((k-1)*MAXLIMAX*MAXLJMAX)			&
	                -(fluxps(k+1) - fluxps(k))*dhs1i(k+1))

	    k = k+1
	  endif

	enddo

      xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX) =                   &
            amax1(0., xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)            &
                      + fluxk(:,KMAX_MID)*dhs1i(KMAX_MID+1))
      ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX) =                         &
            amax1(0., ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)            &
                       + fluxps(KMAX_MID)*dhs1i(KMAX_MID+1))

	end subroutine advvk
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine advvdifvk(xn_adv,ps3d,sdot,skh	&
		,ds3,ds4				&
		,psfac,dt_s)

!     executes advection with a. bott's integreated flux-form
!     using 2'nd order polynomial in the vertical.

	use ModelConstants_ml  , only : KCHEMTOP, EPSIL
	use GenSpec_adv_ml , only : NSPEC_ADV
!hf u2      use My_Runmode_ml , only : ADVEC_TYPE
	implicit none

!	input
      real,intent(in)::  sdot(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1)
      real,intent(in)::  skh(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1)
      real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)
	real,intent(in)::  psfac,dt_s

!	output
      real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)
      real ,intent(inout):: ps3d(0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)

!	local
      real fluxk(NSPEC_ADV,KMAX_MID),fluxps(KMAX_MID),fc(KMAX_MID)

!	local

	integer  k,n,k1
	integer klimlow,klimhig,n1k
	real totk(NSPEC_ADV),totps
	real zzfl1,zzfl2,zzfl3
      real adif(KMAX_MID),bdif(KMAX_MID),cdif(KMAX_MID)            &
            ,e1(KMAX_MID),f1(NSPEC_ADV,KMAX_MID),fps(KMAX_MID)
	real fc1,fc2,fc3

      do k = 1,KMAX_MID-1
	  fc(k) = sdot(k*MAXLIMAX*MAXLJMAX)*dt_s
	  adif(k) = skh(k*MAXLIMAX*MAXLJMAX)*ds3(k)
	  bdif(k+1) = skh(k*MAXLIMAX*MAXLJMAX)*ds4(k)
	enddo
      fc(KMAX_MID) = -1.

!-------------- calculate the advection ----------------------------

	klimlow = 2
	if(fc(1).ge.0.)klimlow=3
      klimhig = KMAX_MID
      if(fc(KMAX_MID-1).lt.0)klimhig = KMAX_MID-1

	fluxk(:,1) = 0.
	fluxps(1) = 0.

	if(fc(1).ge.0.)then

	  fc1 = fc(1)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
	  zzfl2 = alfbegnew(1)*fc1		&
		 + alfbegnew(2)*fc2		&
		 + alfbegnew(3)*fc3
	  zzfl3 = alfnew(7,2,0)*fc1		&
		 + alfnew(8,2,0)*fc2		&
		 + alfnew(9,2,0)*fc3

	  fluxk(:,2) = amax1(0.,xn_adv(:,0)*zzfl2	&
		+xn_adv(:,MAXLIMAX*MAXLJMAX)*zzfl3)
	  fluxps(2) = amax1(0.,ps3d(0)*zzfl2		&
		+ps3d(MAXLIMAX*MAXLJMAX)*zzfl3)

	endif

	do k = klimlow,klimhig
	  fc1 = fc(k-1)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
	  n1k = 0
	  if(fc1.lt.0.)n1k=1

	  zzfl1 = alfnew(1,k,n1k)*fc1		&
		 + alfnew(2,k,n1k)*fc2		&
		 + alfnew(3,k,n1k)*fc3
	  zzfl2 = alfnew(4,k,n1k)*fc1		&
		 + alfnew(5,k,n1k)*fc2		&
		 + alfnew(6,k,n1k)*fc3
	  zzfl3 = alfnew(7,k,n1k)*fc1		&
		 + alfnew(8,k,n1k)*fc2		&
		 + alfnew(9,k,n1k)*fc3

	  k1 = k-2+n1k

	  fluxk(:,k) = amax1(0.,xn_adv(:,(k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1 &
		+xn_adv(:,k1*MAXLIMAX*MAXLJMAX)*zzfl2			 &
		+xn_adv(:,(k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)
	  fluxps(k) = amax1(0.,ps3d((k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1	 &
		+ps3d(k1*MAXLIMAX*MAXLJMAX)*zzfl2			 &
		+ps3d((k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)

	enddo

      if(fc(KMAX_MID-1).lt.0.)then

        fc1 = fc(KMAX_MID-1)
	  fc2 = fc1*fc1
	  fc3 = fc1*fc2
        zzfl1 = alfnew(1,KMAX_MID,1)*fc1            &
             + alfnew(2,KMAX_MID,1)*fc2      &
             + alfnew(3,KMAX_MID,1)*fc3
	  zzfl2 = alfendnew(1)*fc1		&
		 + alfendnew(2)*fc2		&
		 + alfendnew(3)*fc3

        fluxk(:,KMAX_MID) = amax1(0.                        &
            ,xn_adv(:,(KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1      &
            +xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2) 
        fluxps(KMAX_MID) = amax1(0.                        &
            ,ps3d((KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1      &
            +ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2) 

	endif

	k=2

	do while(.true.)

	  do while(fc(k-1).ge.0.)
              fluxk(:,k) = 						&
		amin1(xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX)*dhs1(k),	&
			fluxk(:,k))
	      f1(:,k-1) = amax1(0.,xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX) -  &
			(fluxk(:,k) - fluxk(:,k-1))*dhs1i(k))
              fluxps(k) = 						&
		amin1(ps3d((k-2)*MAXLIMAX*MAXLJMAX)*dhs1(k),		&
			fluxps(k))
	      fps(k-1) = amax1(0.,ps3d((k-2)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxps(k) - fluxps(k-1))*dhs1i(k))
	    k=k+1
          if(k.gt.KMAX_MID)goto 435
	  enddo

	  do while(fc(k).lt.0.)
	      fluxk(:,k) = -amin1(xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)	&
			*dhs1(k+1),fluxk(:,k))
	      f1(:,k-1) = amax1(0.,xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxk(:,k) - fluxk(:,k-1))*dhs1i(k))
	      fluxps(k) = -amin1(ps3d((k-1)*MAXLIMAX*MAXLJMAX)		&
			*dhs1(k+1),fluxps(k))
	      fps(k-1) = amax1(0.,ps3d((k-2)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxps(k) - fluxps(k-1))*dhs1i(k))

	    k=k+1
          if(k.gt.KMAX_MID)goto 435
	  enddo
	  totk(:) = amin1(xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1)/ 	&
		(fluxk(:,k) 						&
		+ fluxk(:,k+1)+ EPSIL),1.)
	  fluxk(:,k) = -fluxk(:,k)*totk(:)
	  fluxk(:,k+1) = fluxk(:,k+1)*totk(:)
	  f1(:,k-1) = amax1(0.,xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxk(:,k) - fluxk(:,k-1))*dhs1i(k))
	  f1(:,k) = amax1(0.,xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))

	  totps = amin1(ps3d((k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1)/	&
		(fluxps(k)+ fluxps(k+1)+ EPSIL),1.)

	  fluxps(k) = -fluxps(k)*totps
	  fluxps(k+1) = fluxps(k+1)*totps
	  fps(k-1) = amax1(0.,ps3d((k-2)*MAXLIMAX*MAXLJMAX) - 	&
			(fluxps(k) - fluxps(k-1))*dhs1i(k))
	  fps(k) = amax1(0.,ps3d((k-1)*MAXLIMAX*MAXLJMAX) - 		&
			(fluxps(k+1) - fluxps(k))*dhs1i(k+1))
	  k = k+2

        if(k.gt.KMAX_MID)goto 435

	enddo

435	continue


	if(ADVEC_TYPE==1)then
          fps(KMAX_MID) =                               &
            amax1(0.,ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)      &
                  + fluxps(KMAX_MID)*dhs1i(KMAX_MID+1))
	endif
	if(ADVEC_TYPE==0)then
	    fps(:) = psfac
	endif

      cdif(KMAX_MID) = 1. + bdif(KMAX_MID)
      e1(KMAX_MID) = bdif(KMAX_MID)/cdif(KMAX_MID)
      f1(:,KMAX_MID) = amax1(0.,xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)      &
                  + fluxk(:,KMAX_MID)*dhs1i(KMAX_MID+1))            &
                   /cdif(KMAX_MID)/fps(KMAX_MID)

      do k = KMAX_MID-1,2,-1
     
	  cdif(k) = 1. + bdif(k) + adif(k) - adif(k)*e1(k+1)
	  e1(k) = bdif(k)/cdif(k)
	  f1(:,k) = (f1(:,k)/fps(k) + adif(k)*f1(:,k+1))	&
			/cdif(k)
    
	enddo

!pw	if(KCHEMTOP.eq.1)then
	  cdif(1) = 1. + adif(1) - adif(1)*e1(2)
          xn_adv(:,0) = (f1(:,1)/fps(1) + adif(1)*f1(:,2))	&
			/cdif(1)
!pw	else
!pw	  cdif(1) = 1. + adif(1) - adif(1)*e1(2)
!pw          f1(:,1) = (f1(:,1)/fps(1) + adif(1)*f1(:,2))	&
!pw			/cdif(1)
!pw          xn_adv(:,MAXLIMAX*MAXLJMAX) = e1(2)*f1(:,1)		&
!pw			+ f1(:,2)
          xn_adv(:,MAXLIMAX*MAXLJMAX) = e1(2)*xn_adv(:,0)	&
			+ f1(:,2)
!pw	endif

!pw      do k = KCHEMTOP+1,KMAX_MID
        do k = 3,KMAX_MID
          xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) = 			&
			e1(k)*xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX)	&
			+ f1(:,k)
	enddo

	end subroutine advvdifvk
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
	subroutine advx(vel,velbeg,velend	&
		,xn_adv,xnbeg,xnend		&
		,ps3d,psbeg,psend		&
		,xm2loc,xmdloc			&
		,dth,fac1)

!     executes advection with a.bott's integrated flux method using
!     4'th order polynomials in the y-direction.
!
!     modified by pw february 2002:  Takes into account the mapping factor 
!     in such a way that a Courant number of one corresponds exactly to "empty" a cell. 
!     (small effects on results: less than 1%)

	use Par_ml   , only : li0,li1,limax
	use GenSpec_adv_ml , only : NSPEC_ADV
        use MassBudget_ml , only : fluxin,fluxout
	implicit none

!	parameter:
!	input
	real,intent(in):: vel(0:MAXLIMAX),velbeg, velend
	real,intent(in):: xnbeg(NSPEC_ADV,3)			&
		,xnend(NSPEC_ADV,3)
	real,intent(in):: psbeg(3)				&
		,psend(3)					&
		,xm2loc(0:MAXLIMAX+1)				&
		,xmdloc(0:MAXLIMAX+1)
	real,intent(in):: dth,fac1

!	input+output
	real ,intent(inout)::xn_adv(NSPEC_ADV,MAXLIMAX)
	real ,intent(inout)::ps3d(MAXLIMAX)

!      output fluxin,fluxout

!	local

	integer ij, ijn,ijll
	integer limtlow,limthig
	integer lijb,lije
	real ijn1
	real x1, x2, x3,hh3,hh4
	real y0,y1,y2,y3
	real zzfc(5,-1:MAXLIMAX+1)
	real fc(-1:MAXLIMAX+1)
	real flux(NSPEC_ADV,-1:MAXLIMAX+1)
	real fluxps(-1:MAXLIMAX+1)
	real hel1(NSPEC_ADV),hel2(NSPEC_ADV)
	real hel1ps,hel2ps
	integer ijpasses
	integer ijb1(MAXLIMAX),ije1(MAXLIMAX)
	integer ijb2(MAXLIMAX),ije2(MAXLIMAX),ijb3(MAXLIMAX)
	logical ijdoend

!-----------------------------------------------------------------------

	limtlow = li0-1
	if (li0.eq.1) then
	  if (vel(0) .gt. 0..and.velbeg.lt.0.) then
	    fc(-1) = velbeg*dth
	    limtlow = -1

	    y0 = fc(-1)
	    x1 = 1.+2.*y0*xm2loc(0)
	    x2 = x1*x1
	    y3 = xmdloc(0)*(1.-x2)/3840.
	    y1 = 5.*y3
	    y2 = x1*y3
	    hh3 = (116.-4.*x2)*y2
	    hh4 = (2.*x2-66.)*y1
	    zzfc(3,-1) = - y0 - (214. - 6.*x2)*y2
	    zzfc(5,-1) = (y2-y1)*(x2-9.)
	    zzfc(1,-1) = (y2+y1)*(x2-9.)
	    zzfc(4,-1) = hh3+hh4
	    zzfc(2,-1) = hh3-hh4

	  endif
	endif

	do 10 ij = li0-1,li1
	  fc(ij) = vel(ij)*dth

	  ijn1 = sign(1.,fc(ij))
          ijn = ij + nint(0.5*(1-ijn1))
          
	  y0 = ijn1*fc(ij)
	  x1 = 1.-2.*y0*xm2loc(ijn)
	  x2 = x1*x1
	  y3 = xmdloc(ijn)*(1.-x2)/3840.
	  y1 = 5.*ijn1*y3
	  y2 = x1*y3
	  hh3 = (116.-4.*x2)*y2
	  hh4 = (66.-2.*x2)*y1
	  zzfc(3,ij) = y0 - (214. - 6.*x2)*y2
	  zzfc(5,ij) = (y2+y1)*(x2-9.)
	  zzfc(1,ij) = (y2-y1)*(x2-9.)
	  zzfc(4,ij) = hh3+hh4
	  zzfc(2,ij) = hh3-hh4

10	continue

	limthig = li1
	if (li1.eq.limax) then
	  if (vel(li1).lt.0..and.velend.gt.0.)then
	    fc(li1+1) = velend*dth
	    limthig = li1+1

	    y0 = fc(li1+1)
	    x1 = 1.-2.*y0*xm2loc(li1+1)
	    x2 = x1*x1
	    y3 = xmdloc(li1+1)*(1.-x2)/3840.
	    y1 = 5.*y3
	    y2 = x1*y3
	    hh3 = (116.-4.*x2)*y2
	    hh4 = (66.-2.*x2)*y1
	    zzfc(3,li1+1) = y0 - (214.-6.*x2)*y2
	    zzfc(5,li1+1) = (y2+y1)*(x2-9.)
	    zzfc(1,li1+1) = (y2-y1)*(x2-9.)
	    zzfc(4,li1+1) = hh3+hh4
	    zzfc(2,li1+1) = hh3-hh4

	  endif
	endif

!------- boundary treatment -----------------------------------------   

!        helping values at the boundaries are found by linear
!        extrapolation in cases of outflow, and by assuming constant 
!        values in inflow cases.

!        calculate the coefficients in the polynomial, the
!        normalized fluxes, and limit them for positivness

	if(limtlow.eq.-1)then

!     integrated flux form

	  flux(:,-1) = max(0.,xn_adv(:,2)*zzfc(5,-1)		&
		+ xn_adv(:,1)*zzfc(4,-1)				&
		+ xnbeg(:,3)*zzfc(3,-1)				&
		+ xnbeg(:,2)*zzfc(2,-1) + xnbeg(:,1)*zzfc(1,-1))
	  flux(:,0) = max(0.,xn_adv(:,2)*zzfc(5,0)		&
		+ xn_adv(:,1)*zzfc(4,0)				&
		+ xnbeg(:,3)*zzfc(3,0)				&
		+ xnbeg(:,2)*zzfc(2,0) + xnbeg(:,1)*zzfc(1,0))
	  fluxps(-1) = max(0.,ps3d(2)*zzfc(5,-1)			&
		+ ps3d(1)*zzfc(4,-1)				&
		+ psbeg(3)*zzfc(3,-1)				&
		+ psbeg(2)*zzfc(2,-1) + psbeg(1)*zzfc(1,-1))
	  fluxps(0) = max(0.,ps3d(2)*zzfc(5,0)			&
		+ ps3d(1)*zzfc(4,0)				&
		+ psbeg(3)*zzfc(3,0)				&
		+ psbeg(2)*zzfc(2,0) + psbeg(1)*zzfc(1,0))

	else

!     integrated flux form

	  if(fc(li0-1).ge.0.)then
	    flux(:,li0-1) = max(0.,xn_adv(:,li0+1)*zzfc(5,li0-1)		&
		+ xn_adv(:,li0)*zzfc(4,li0-1)				&
		+ xnbeg(:,3)*zzfc(3,li0-1)				&
		+ xnbeg(:,2)*zzfc(2,li0-1) + xnbeg(:,1)*zzfc(1,li0-1))
	    fluxps(li0-1) = max(0.,ps3d(li0+1)*zzfc(5,li0-1)		&
		+ ps3d(li0)*zzfc(4,li0-1)					&
		+ psbeg(3)*zzfc(3,li0-1)					&
		+ psbeg(2)*zzfc(2,li0-1) + psbeg(1)*zzfc(1,li0-1))
	  else
	    flux(:,li0-1) = max(0.,xn_adv(:,li0+2)*zzfc(5,li0-1)		&
		+ xn_adv(:,li0+1)*zzfc(4,li0-1)				&
		+ xn_adv(:,li0)*zzfc(3,li0-1)				&
		+ xnbeg(:,3)*zzfc(2,li0-1) + xnbeg(:,2)*zzfc(1,li0-1))
	    fluxps(li0-1) = max(0.,ps3d(li0+2)*zzfc(5,li0-1)		&
		+ ps3d(li0+1)*zzfc(4,li0-1)				&
		+ ps3d(li0)*zzfc(3,li0-1)					&
		+ psbeg(3)*zzfc(2,li0-1) + psbeg(2)*zzfc(1,li0-1))
	  endif
	endif

!     integrated flux form

	if(fc(li0).ge.0.)then
	  flux(:,li0) = max(0.,xn_adv(:,li0+2)*zzfc(5,li0)		&
		+ xn_adv(:,li0+1)*zzfc(4,li0)				&
		+ xn_adv(:,li0)*zzfc(3,li0)				&
		+ xnbeg(:,3)*zzfc(2,li0) + xnbeg(:,2)*zzfc(1,li0))
	  fluxps(li0) = max(0.,ps3d(li0+2)*zzfc(5,li0)			&
		+ ps3d(li0+1)*zzfc(4,li0)					&
		+ ps3d(li0)*zzfc(3,li0)					&
		+ psbeg(3)*zzfc(2,li0) + psbeg(2)*zzfc(1,li0))
	else
	  flux(:,li0) = max(0.,xn_adv(:,li0+3)*zzfc(5,li0)		&
		+ xn_adv(:,li0+2)*zzfc(4,li0)				&
		+ xn_adv(:,li0+1)*zzfc(3,li0)				&
		+ xn_adv(:,li0)*zzfc(2,li0) + xnbeg(:,3)*zzfc(1,li0))
	  fluxps(li0) = max(0.,ps3d(li0+3)*zzfc(5,li0)			&
		+ ps3d(li0+2)*zzfc(4,li0)					&
		+ ps3d(li0+1)*zzfc(3,li0)					&
		+ ps3d(li0)*zzfc(2,li0) + psbeg(3)*zzfc(1,li0))
	endif

	if(fc(li0+1).ge.0.)then

!     integrated flux form

	  flux(:,li0+1) = max(0.,xn_adv(:,li0+3)*zzfc(5,li0+1)		&
		+ xn_adv(:,li0+2)*zzfc(4,li0+1)				&
		+ xn_adv(:,li0+1)*zzfc(3,li0+1)				&
		+ xn_adv(:,li0)*zzfc(2,li0+1) + xnbeg(:,3)*zzfc(1,li0+1))
	  fluxps(li0+1) = max(0.,ps3d(li0+3)*zzfc(5,li0+1)		&
		+ ps3d(li0+2)*zzfc(4,li0+1)				&
		+ ps3d(li0+1)*zzfc(3,li0+1)				&
		+ ps3d(li0)*zzfc(2,li0+1) + psbeg(3)*zzfc(1,li0+1))
	endif

	lijb = li0+2
	if(fc(li0+1).lt.0.)lijb = li0+1
	lije = li1-3
	if(fc(li1-2).ge.0.)lije = li1-2

	do ij = lijb,lije

	  ijn1 = sign(1.,fc(ij))

!     integrated flux form

	  ijn = ij+nint(0.5*(1.-ijn1))
	  flux(:,ij) = max(0.,xn_adv(:,ijn+2)*zzfc(5,ij)		&
		+ xn_adv(:,ijn+1)*zzfc(4,ij)			&
		+ xn_adv(:,ijn)*zzfc(3,ij)			&
		+ xn_adv(:,ijn-1)*zzfc(2,ij)			& 
		+ xn_adv(:,ijn-2)*zzfc(1,ij))
	  fluxps(ij) = max(0.,ps3d(ijn+2)*zzfc(5,ij)	&
		+ ps3d(ijn+1)*zzfc(4,ij)			&
		+ ps3d(ijn)*zzfc(3,ij)			&
		+ ps3d(ijn-1)*zzfc(2,ij)			& 
		+ ps3d(ijn-2)*zzfc(1,ij))

	enddo

	if(fc(li1-2).lt.0)then

!     integrated flux form

	  flux(:,li1-2) = max(0.,xnend(:,1)*zzfc(5,li1-2)	&
		+ xn_adv(:,li1)*zzfc(4,li1-2)			&
		+ xn_adv(:,li1-1)*zzfc(3,li1-2)			&
		+ xn_adv(:,li1-2)*zzfc(2,li1-2)			& 
		+ xn_adv(:,li1-3)*zzfc(1,li1-2))
	  fluxps(li1-2) = max(0.,psend(1)*zzfc(5,li1-2)		&
		+ ps3d(li1)*zzfc(4,li1-2)				&
		+ ps3d(li1-1)*zzfc(3,li1-2)			&
		+ ps3d(li1-2)*zzfc(2,li1-2)			& 
		+ ps3d(li1-3)*zzfc(1,li1-2))
	endif

!     integrated flux form

	if(fc(li1-1).ge.0)then
	  flux(:,li1-1) = max(0.,xnend(:,1)*zzfc(5,li1-1)	&
		+ xn_adv(:,li1)*zzfc(4,li1-1)			&
		+ xn_adv(:,li1-1)*zzfc(3,li1-1)			&
		+ xn_adv(:,li1-2)*zzfc(2,li1-1)			& 
		+ xn_adv(:,li1-3)*zzfc(1,li1-1))
	  fluxps(li1-1) = max(0.,psend(1)*zzfc(5,li1-1)	&
		+ ps3d(li1)*zzfc(4,li1-1)				&
		+ ps3d(li1-1)*zzfc(3,li1-1)			&
		+ ps3d(li1-2)*zzfc(2,li1-1)			& 
		+ ps3d(li1-3)*zzfc(1,li1-1))
	else
	  flux(:,li1-1) = max(0.,xnend(:,2)*zzfc(5,li1-1)	&
		+ xnend(:,1)*zzfc(4,li1-1)			&
		+ xn_adv(:,li1)*zzfc(3,li1-1)			&
		+ xn_adv(:,li1-1)*zzfc(2,li1-1)			& 
		+ xn_adv(:,li1-2)*zzfc(1,li1-1))
	  fluxps(li1-1) = max(0.,psend(2)*zzfc(5,li1-1)		&
		+ psend(1)*zzfc(4,li1-1)				&
		+ ps3d(li1)*zzfc(3,li1-1)				&
		+ ps3d(li1-1)*zzfc(2,li1-1)			& 
		+ ps3d(li1-2)*zzfc(1,li1-1))
	endif

!     integrated flux form
	if(limthig.eq.li1)then

	  if(fc(li1).ge.0)then
	    flux(:,li1) = max(0.,xnend(:,2)*zzfc(5,li1)		&
		+ xnend(:,1)*zzfc(4,li1)				&
		+ xn_adv(:,li1)*zzfc(3,li1)			&
		+ xn_adv(:,li1-1)*zzfc(2,li1)			& 
		+ xn_adv(:,li1-2)*zzfc(1,li1))
	    fluxps(li1) = max(0.,psend(2)*zzfc(5,li1)		&
		+ psend(1)*zzfc(4,li1)				&
		+ ps3d(li1)*zzfc(3,li1)				&
		+ ps3d(li1-1)*zzfc(2,li1)				& 
		+ ps3d(li1-2)*zzfc(1,li1))
	  else
	    flux(:,li1) = max(0.,xnend(:,3)*zzfc(5,li1)		&
		+ xnend(:,2)*zzfc(4,li1)				&
		+ xnend(:,1)*zzfc(3,li1)				&
		+ xn_adv(:,li1)*zzfc(2,li1)			& 
		+ xn_adv(:,li1-1)*zzfc(1,li1))
	    fluxps(li1) = max(0.,psend(3)*zzfc(5,li1)		&
		+ psend(2)*zzfc(4,li1)				&
		+ psend(1)*zzfc(3,li1)				&
		+ ps3d(li1)*zzfc(2,li1)				& 
		+ ps3d(li1-1)*zzfc(1,li1))
	  endif

	else

!     integrated flux form

	  flux(:,li1) = max(0.,xnend(:,3)*zzfc(5,li1)		&
		+ xnend(:,2)*zzfc(4,li1)				&
		+ xnend(:,1)*zzfc(3,li1)				&
		+ xn_adv(:,li1)*zzfc(2,li1)			& 
		+ xn_adv(:,li1-1)*zzfc(1,li1))
	  flux(:,li1+1) = max(0.,xnend(:,3)*zzfc(5,li1+1)		&
		+ xnend(:,2)*zzfc(4,li1+1)			&
		+ xnend(:,1)*zzfc(3,li1+1)			&
		+ xn_adv(:,li1)*zzfc(2,li1+1)			& 
		+ xn_adv(:,li1-1)*zzfc(1,li1+1))
	  fluxps(li1) = max(0.,psend(3)*zzfc(5,li1)		&
		+ psend(2)*zzfc(4,li1)				&
		+ psend(1)*zzfc(3,li1)				&
		+ ps3d(li1)*zzfc(2,li1)				& 
		+ ps3d(li1-1)*zzfc(1,li1))
	  fluxps(li1+1) = max(0.,psend(3)*zzfc(5,li1+1)		&
		+ psend(2)*zzfc(4,li1+1)				&
		+ psend(1)*zzfc(3,li1+1)				&
		+ ps3d(li1)*zzfc(2,li1+1)				& 
		+ ps3d(li1-1)*zzfc(1,li1+1))

	endif


	if(limtlow.eq.-1)then
		hel1(:) = xnbeg(:,3)*xmdloc(0)
		hel2(:) = flux(:,0) +  flux(:,-1)
		where(hel1(:).lt.hel2(:))flux(:,0)		&
			= flux(:,0)*hel1(:)/hel2(:)
		hel1ps = psbeg(3)*xmdloc(0)
		hel2ps = fluxps(0) +  fluxps(-1)
		if(hel1ps.lt.hel2ps)fluxps(0) = fluxps(0)*hel1ps/hel2ps
	    	ij = 1
	else
	  if(fc(li0-1).ge.0.) then
	    flux(:,li0-1) = amin1(xnbeg(:,3)		&
			*xmdloc(li0-1)			&
			,flux(:,li0-1))
	    fluxps(li0-1) = amin1(psbeg(3)		&
			*xmdloc(li0-1)			&
			,fluxps(li0-1))
	    ij = li0
	  else
	    if(fc(li0).lt.0.) then
	      flux(:,li0-1) = -amin1(xn_adv(:,li0)	&
			*xmdloc(li0),flux(:,li0-1))
	      fluxps(li0-1) = -amin1(ps3d(li0)		&
			*xmdloc(li0),fluxps(li0-1))
	      ij = li0
	    else
		hel1(:) = xn_adv(:,li0)*xmdloc(li0)
		hel2(:) = flux(:,li0) + flux(:,li0-1)
		where(hel1(:).lt.hel2(:))
		  flux(:,li0-1) = -flux(:,li0-1)*hel1(:)/hel2(:)
		  flux(:,li0) = flux(:,li0)*hel1(:)/hel2(:)
		  xn_adv(:,li0) = 0.
		else where
		  flux(:,li0-1) = -flux(:,li0-1)
		  xn_adv(:,li0) =xm2loc(li0)*(hel1(:)-hel2(:))
		end where
		hel1ps = ps3d(li0)*xmdloc(li0)
		hel2ps = fluxps(li0) +  fluxps(li0-1)
		if(hel1ps.lt.hel2ps)then
		  fluxps(li0-1) = - fluxps(li0-1)*hel1ps/hel2ps
		  fluxps(li0) = fluxps(li0)*hel1ps/hel2ps
		  ps3d(li0) = 0.
		else
		  fluxps(li0-1) = -fluxps(li0-1)
		  ps3d(li0) =xm2loc(li0)*(hel1ps-hel2ps)
		endif
	    	ij = li0+1
	    endif
	  endif
	endif

	ijpasses = 0
	do while(.true.)

	  ijpasses = ijpasses+1
	  ijb1(ijpasses) = ij
	  ije1(ijpasses) = -5
	  do while(fc(ij).ge.0.)
	    ije1(ijpasses) = ij
	    ij = ij+1
	    if(ij.gt.li1-1)then
	      ijb2(ijpasses) = ij
	      ije2(ijpasses) = -5
	      ijb3(ijpasses) = -5
	      goto 257
	    endif
	  enddo
	  ijb2(ijpasses) = ij
	  ije2(ijpasses) = -5
	  do while(fc(ij+1).lt.0.)
	    ije2(ijpasses) = ij
	    ij = ij+1
	    if(ij.gt.li1-1)then
	      ijb3(ijpasses) = -5
	      goto 257
	    endif
	  enddo
	  ijb3(ijpasses) = ij
	  ij = ij+2
	  if(ij.gt.li1-1)goto 257
	enddo

257	continue
	ijdoend = .false.
	if(ij.eq.li1)ijdoend=.true.

	do ijll = 1,ijpasses

	  do ij = ijb1(ijll),ije1(ijll)
	      flux(:,ij) = amin1(xn_adv(:,ij)*xmdloc(ij)	&
			,flux(:,ij))
	      xn_adv(:,ij) =					&
			amax1(0.,xn_adv(:,ij)			&
			-xm2loc(ij)				&
			*(flux(:,ij) - flux(:,ij-1)))
	      fluxps(ij) = amin1(ps3d(ij)*xmdloc(ij)		&
			,fluxps(ij))
	      ps3d(ij) =					&
			amax1(0.,ps3d(ij)			&
			-xm2loc(ij)				&
			*(fluxps(ij) - fluxps(ij-1)))
	  enddo
	  do ij = ijb2(ijll),ije2(ijll)
	      flux(:,ij) = -amin1(xn_adv(:,ij+1)*xmdloc(ij+1)	&
			,flux(:,ij))
	      xn_adv(:,ij) =					&
			amax1(0.,xn_adv(:,ij)			&
			-xm2loc(ij)*(flux(:,ij)			&
			- flux(:,ij-1)))
	      fluxps(ij) = -amin1(ps3d(ij+1)*xmdloc(ij+1)	&
			,fluxps(ij))
	      ps3d(ij) =					&
			amax1(0.,ps3d(ij)			&
			-xm2loc(ij)*(fluxps(ij)			&
			- fluxps(ij-1)))
	  enddo
	  ij = ijb3(ijll)
	  if(ij.lt.-3) goto 357
	    hel1(:) = xn_adv(:,ij+1)*xmdloc(ij+1)
	    hel2(:) = flux(:,ij+1) +  flux(:,ij)
	    where(hel1(:).lt.hel2(:))
	      flux(:,ij) = - flux(:,ij)*hel1(:)/hel2(:)
	      flux(:,ij+1) = flux(:,ij+1)*hel1(:)/hel2(:)
	      xn_adv(:,ij+1) = 0.
	    else where
	      flux(:,ij) = -flux(:,ij)
	      xn_adv(:,ij+1) = xm2loc(ij+1)*(hel1(:)-hel2(:))
	    end where
	    xn_adv(:,ij) =					&
			amax1(0.,xn_adv(:,ij)			&
			-xm2loc(ij)*(flux(:,ij) - flux(:,ij-1)))
	    hel1ps = ps3d(ij+1)*xmdloc(ij+1)
	    hel2ps = fluxps(ij+1) +  fluxps(ij)
	    if(hel1ps.lt.hel2ps)then
	      fluxps(ij) = -fluxps(ij)*hel1ps/hel2ps
	      fluxps(ij+1) = fluxps(ij+1)*hel1ps/hel2ps
	      ps3d(ij+1) = 0.
	    else
	      fluxps(ij) = -fluxps(ij)
	      ps3d(ij+1) = xm2loc(ij+1)*(hel1ps-hel2ps)
	    endif
	    ps3d(ij) =						&
			amax1(0.,ps3d(ij)			&
			-xm2loc(ij)*(fluxps(ij) - fluxps(ij-1)))
	enddo

357	continue

	if(ijdoend)then
	  if(limthig.eq.li1+1)then

	    hel1(:) = xnend(:,1)*xmdloc(li1+1)
	    hel2(:) = flux(:,li1+1) + flux(:,li1)
	    where(hel1(:).lt.hel2(:))
	      flux(:,li1) = 					&
			- flux(:,li1)*hel1(:)/hel2(:)
	    else where
	      flux(:,li1) = -flux(:,li1)
	    end where
	    xn_adv(:,li1) =amax1(0.				&
			,xn_adv(:,li1)-xm2loc(li1)		&
			*(flux(:,li1)- flux(:,li1-1)))
	    hel1ps = psend(1)*xmdloc(li1+1)
	    hel2ps = fluxps(li1+1) + fluxps(li1)
	    if(hel1ps.lt.hel2ps)then
	      fluxps(li1) = 					&
			-fluxps(li1)*hel1ps/hel2ps
	    else
	      fluxps(li1) = -fluxps(li1)
	    endif
	    ps3d(li1) =amax1(0.					&
			,ps3d(li1)-xm2loc(li1)			&
			*(fluxps(li1)- fluxps(li1-1)))

	  else

	    if(fc(li1).ge.0.) then
	      flux(:,li1) = amin1(xn_adv(:,li1)			&
			*xmdloc(li1),flux(:,li1))
	      xn_adv(:,li1) =amax1(0.				&
			,xn_adv(:,li1)				&
			-xm2loc(li1)*(flux(:,li1)- flux(:,li1-1)))
	      fluxps(li1) = amin1(ps3d(li1)			&
			*xmdloc(li1),fluxps(li1))
	      ps3d(li1) =amax1(0.				&
			,ps3d(li1)				&
			-xm2loc(li1)*(fluxps(li1)- fluxps(li1-1)))
	    else
		flux(:,li1) = -amin1(xnend(:,1)*xmdloc(li1+1)	&
			,flux(:,li1))
		xn_adv(:,li1) =amax1(0.				&
			,xn_adv(:,li1)-xm2loc(li1)		&
			*(flux(:,li1)- flux(:,li1-1)))
		fluxps(li1) = -amin1(psend(1)*xmdloc(li1+1)	&
			,fluxps(li1))
		ps3d(li1) =amax1(0.				&
			,ps3d(li1)-xm2loc(li1)			&
			*(fluxps(li1)- fluxps(li1-1)))
	    endif
	  endif
	endif

!     accumulation of the boundary fluxes

	if (li0.eq.2) then
	  if(fc(1).ge.0.)then
            fluxin(:) = fluxin(:) + flux(:,1)*fac1
	  else
            fluxout(:) = fluxout(:) - flux(:,1)*fac1
	  endif
	endif

	if (li1.eq.limax-1) then
	  if(fc(li1).ge.0.)then
            fluxout(:) = fluxout(:) + flux(:,li1)*fac1
	  else
            fluxin(:) = fluxin(:) - flux(:,li1)*fac1
	  endif
	endif

	end subroutine advx

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
	subroutine advy(vel,velbeg,velend		&
		,xn_adv,xnbeg,xnend			&
		,ps3d,psbeg,psend			&
		,xm2loc,xmdloc				&
		,dth,fac1)

!     executes advection with a.bott's integrated flux method using
!     4'th order polynomials in the y-direction.
!
!     modified by pw february 2002:  Takes into account the mapping factor 
!     in such a way that a Courant number of one corresponds exactly to "empty" a cell. 
!     (small effects on results: less than 1%)
!

	use Par_ml   , only : lj0,lj1,ljmax
	use GenSpec_adv_ml , only : NSPEC_ADV
        use MassBudget_ml , only : fluxin,fluxout
	implicit none

!	parameter:
!	input
	real,intent(in):: vel(0:MAXLIMAX*MAXLJMAX),velbeg, velend
	real,intent(in):: xnbeg(NSPEC_ADV,3)		&
		,xnend(NSPEC_ADV,3)
	real,intent(in):: psbeg(3)			&
		,psend(3)				&
		,xm2loc(0:MAXLJMAX+1)			&
		,xmdloc(0:MAXLJMAX+1)
	real,intent(in):: dth,fac1

!	input+output
	real ,intent(inout)::xn_adv(NSPEC_ADV,MAXLIMAX:MAXLIMAX*MAXLJMAX)
	real ,intent(inout)::ps3d(MAXLIMAX:MAXLIMAX*MAXLJMAX)

!      output fluxin,fluxout

!	local

	integer ij, ijn,ijll
	integer limtlow,limthig
	integer lijb,lije
	real ijn1
	real x1, x2, x3,hh3,hh4
	real y0,y1,y2,y3
	real zzfc(5,-1:MAXLJMAX+1)
	real fc(-1:MAXLJMAX+1)
	real flux(NSPEC_ADV,-1:MAXLJMAX+1)
	real fluxps(-1:MAXLJMAX+1)
	real hel1(NSPEC_ADV),hel2(NSPEC_ADV)
	real hel1ps,hel2ps
	integer ijpasses
	integer ijb1(MAXLJMAX),ije1(MAXLJMAX)
	integer ijb2(MAXLJMAX),ije2(MAXLJMAX),ijb3(MAXLJMAX)
	logical ijdoend

!-----------------------------------------------------------------------

	limtlow = lj0-1
	if (lj0.eq.1) then
	  if (vel(0) .gt. 0..and.velbeg.lt.0.) then
	    fc(-1) = velbeg*dth
	    limtlow = -1

	    y0 = fc(-1)
	    x1 = 1.+2.*y0*xm2loc(0)
	    x2 = x1*x1
	    y3 = xmdloc(0)*(1.-x2)/3840.
	    y1 = 5.*y3
	    y2 = x1*y3
	    hh3 = (116.-4.*x2)*y2
	    hh4 = (2.*x2-66.)*y1
	    zzfc(3,-1) = - y0 - (214.-6.*x2)*y2
	    zzfc(5,-1) = (y2-y1)*(x2-9.)
	    zzfc(1,-1) = (y2+y1)*(x2-9.)
	    zzfc(4,-1) = hh3+hh4
	    zzfc(2,-1) = hh3-hh4

	  endif
	endif

	do 10 ij = lj0-1,lj1
	  fc(ij) = vel(ij*MAXLIMAX)*dth

	  ijn1 = sign(1.,fc(ij))
          ijn = ij + nint(0.5*(1-ijn1))

	  y0 = ijn1*fc(ij)
	  x1 = 1.-2.*y0*xm2loc(ijn)
	  x2 = x1*x1
	  y3 = xmdloc(ijn)*(1.-x2)/3840.
	  y1 = 5.*ijn1*y3
	  y2 = x1*y3
	  hh3 = (116.-4.*x2)*y2
	  hh4 = (66.-2.*x2)*y1
	  zzfc(3,ij) = y0 - (214.-6.*x2)*y2
	  zzfc(5,ij) = (y2+y1)*(x2-9.)
	  zzfc(1,ij) = (y2-y1)*(x2-9.)
	  zzfc(4,ij) = hh3+hh4
	  zzfc(2,ij) = hh3-hh4

10	continue

	limthig = lj1
	if (lj1.eq.ljmax) then
	  if (vel(lj1*MAXLIMAX).lt.0..and.velend.gt.0.)then
	    fc(lj1+1) = velend*dth
	    limthig = lj1+1

	    y0 = fc(lj1+1)
	    x1 = 1.-2.*y0*xm2loc(lj1+1)
	    x2 = x1*x1
	    y3 = xmdloc(lj1+1)*(1.-x2)/3840.
	    y1 = 5.*y3
	    y2 = x1*y3
	    hh3 = (116.-4.*x2)*y2
	    hh4 = (66.-2.*x2)*y1
	    zzfc(3,lj1+1) = y0 - (214.-6.*x2)*y2
	    zzfc(5,lj1+1) = (y2+y1)*(x2-9.)
	    zzfc(1,lj1+1) = (y2-y1)*(x2-9.)
	    zzfc(4,lj1+1) = hh3+hh4
	    zzfc(2,lj1+1) = hh3-hh4

	  endif
	endif

!------- boundary treatment -----------------------------------------   

!        helping values at the boundaries are found by linear
!        extrapolation in cases of outflow, and by assuming constant 
!        values in inflow cases.

!        calculate the coefficients in the polynomial, the
!        normalized fluxes, and limit them for positivness

	if(limtlow.eq.-1) then


!     integrated flux form

	  flux(:,-1) = max(0.,xn_adv(:,2*MAXLIMAX)*zzfc(5,-1)	&
		+ xn_adv(:,MAXLIMAX)*zzfc(4,-1)			&
		+ xnbeg(:,3)*zzfc(3,-1)				&
		+ xnbeg(:,2)*zzfc(2,-1) + xnbeg(:,1)*zzfc(1,-1))
	  flux(:,0) = max(0.,xn_adv(:,(2)*MAXLIMAX)*zzfc(5,0)	&
		+ xn_adv(:,1*MAXLIMAX)*zzfc(4,0)			&
		+ xnbeg(:,3)*zzfc(3,0)				&
		+ xnbeg(:,2)*zzfc(2,0) + xnbeg(:,1)*zzfc(1,0))
	  fluxps(-1) = max(0.,ps3d(2*MAXLIMAX)*zzfc(5,-1)		&
		+ ps3d(MAXLIMAX)*zzfc(4,-1)			&
		+ psbeg(3)*zzfc(3,-1)				&
		+ psbeg(2)*zzfc(2,-1) + psbeg(1)*zzfc(1,-1))
	  fluxps(0) = max(0.,ps3d((2)*MAXLIMAX)*zzfc(5,0)		&
		+ ps3d(1*MAXLIMAX)*zzfc(4,0)			&
		+ psbeg(3)*zzfc(3,0)				&
		+ psbeg(2)*zzfc(2,0) + psbeg(1)*zzfc(1,0))

	else
!     integrated flux form

	  if(fc(lj0-1).ge.0.)then
	    flux(:,lj0-1) = max(0.,xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(5,lj0-1)&
		+ xn_adv(:,lj0*MAXLIMAX)*zzfc(4,lj0-1)			&
		+ xnbeg(:,3)*zzfc(3,lj0-1)				&
		+ xnbeg(:,2)*zzfc(2,lj0-1) + xnbeg(:,1)*zzfc(1,lj0-1))
	    fluxps(lj0-1) = max(0.,ps3d((lj0+1)*MAXLIMAX)*zzfc(5,lj0-1)	&
		+ ps3d(lj0*MAXLIMAX)*zzfc(4,lj0-1)			&
		+ psbeg(3)*zzfc(3,lj0-1)					&
		+ psbeg(2)*zzfc(2,lj0-1) + psbeg(1)*zzfc(1,lj0-1))
	  else
	    flux(:,lj0-1) = max(0.,xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(5,lj0-1)&
		+ xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(4,lj0-1)		&
		+ xn_adv(:,lj0*MAXLIMAX)*zzfc(3,lj0-1)			&
		+ xnbeg(:,3)*zzfc(2,lj0-1) + xnbeg(:,2)*zzfc(1,lj0-1))
	    fluxps(lj0-1) = max(0.,ps3d((lj0+2)*MAXLIMAX)*zzfc(5,lj0-1)	&
		+ ps3d((lj0+1)*MAXLIMAX)*zzfc(4,lj0-1)			&
		+ ps3d(lj0*MAXLIMAX)*zzfc(3,lj0-1)			&
		+ psbeg(3)*zzfc(2,lj0-1) + psbeg(2)*zzfc(1,lj0-1))
	  endif
	endif

!     integrated flux form

	if(fc(lj0).ge.0.)then
	  flux(:,lj0) = max(0.,xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(5,lj0)	&
		+ xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(4,lj0)			&
		+ xn_adv(:,lj0*MAXLIMAX)*zzfc(3,lj0)			&
		+ xnbeg(:,3)*zzfc(2,lj0) + xnbeg(:,2)*zzfc(1,lj0))
	  fluxps(lj0) = max(0.,ps3d((lj0+2)*MAXLIMAX)*zzfc(5,lj0)	&
		+ ps3d((lj0+1)*MAXLIMAX)*zzfc(4,lj0)			&
		+ ps3d(lj0*MAXLIMAX)*zzfc(3,lj0)				&
		+ psbeg(3)*zzfc(2,lj0) + psbeg(2)*zzfc(1,lj0))
	else
	  flux(:,lj0) = max(0.,xn_adv(:,(lj0+3)*MAXLIMAX)*zzfc(5,lj0)	&
		+ xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(4,lj0)			&
		+ xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(3,lj0)			&
		+ xn_adv(:,lj0*MAXLIMAX)*zzfc(2,lj0) + xnbeg(:,3)*zzfc(1,lj0))
	  fluxps(lj0) = max(0.,ps3d((lj0+3)*MAXLIMAX)*zzfc(5,lj0)	&
		+ ps3d((lj0+2)*MAXLIMAX)*zzfc(4,lj0)			&
		+ ps3d((lj0+1)*MAXLIMAX)*zzfc(3,lj0)			&
		+ ps3d(lj0*MAXLIMAX)*zzfc(2,lj0) + psbeg(3)*zzfc(1,lj0))
	endif

	if(fc(lj0+1).ge.0.)then

!     integrated flux form

	  flux(:,lj0+1) = max(0.,xn_adv(:,(lj0+3)*MAXLIMAX)*zzfc(5,lj0+1)&
		+ xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(4,lj0+1)		&
		+ xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(3,lj0+1)		&
		+ xn_adv(:,lj0*MAXLIMAX)*zzfc(2,lj0+1) + xnbeg(:,3)*zzfc(1,lj0+1))
	  fluxps(lj0+1) = max(0.,ps3d((lj0+3)*MAXLIMAX)*zzfc(5,lj0+1)	&
		+ ps3d((lj0+2)*MAXLIMAX)*zzfc(4,lj0+1)			&
		+ ps3d((lj0+1)*MAXLIMAX)*zzfc(3,lj0+1)			&
		+ ps3d(lj0*MAXLIMAX)*zzfc(2,lj0+1) + psbeg(3)*zzfc(1,lj0+1))
	endif


	lijb = lj0+2
	if(fc(lj0+1).lt.0.)lijb = lj0+1
	lije = lj1-3
	if(fc(lj1-2).ge.0.)lije = lj1-2

	do ij = lijb,lije

	  ijn1 = sign(1.,fc(ij))

!     integrated flux form

	  ijn = ij+nint(0.5*(1.-ijn1))
	  flux(:,ij) = max(0.,xn_adv(:,(ijn+2)*MAXLIMAX)*zzfc(5,ij)	&
		+ xn_adv(:,(ijn+1)*MAXLIMAX)*zzfc(4,ij)			&
		+ xn_adv(:,ijn*MAXLIMAX)*zzfc(3,ij)			&
		+ xn_adv(:,(ijn-1)*MAXLIMAX)*zzfc(2,ij)			& 
		+ xn_adv(:,(ijn-2)*MAXLIMAX)*zzfc(1,ij))
	  fluxps(ij) = max(0.,ps3d((ijn+2)*MAXLIMAX)*zzfc(5,ij)		&
		+ ps3d((ijn+1)*MAXLIMAX)*zzfc(4,ij)			&
		+ ps3d(ijn*MAXLIMAX)*zzfc(3,ij)				&
		+ ps3d((ijn-1)*MAXLIMAX)*zzfc(2,ij)			& 
		+ ps3d((ijn-2)*MAXLIMAX)*zzfc(1,ij))

	enddo

	if(fc(lj1-2).lt.0.)then

!     integrated flux form


	  flux(:,lj1-2) = max(0.,xnend(:,1)*zzfc(5,lj1-2)		&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(4,lj1-2)		&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(3,lj1-2)	&
		+ xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(2,lj1-2) 	&
		+ xn_adv(:,(lj1-3)*MAXLIMAX)*zzfc(1,lj1-2))
	  fluxps(lj1-2) = max(0.,psend(1)*zzfc(5,lj1-2)		&
		+ ps3d(lj1*MAXLIMAX)*zzfc(4,lj1-2)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(3,lj1-2)		&
		+ ps3d((lj1-2)*MAXLIMAX)*zzfc(2,lj1-2)		&
		+ ps3d((lj1-3)*MAXLIMAX)*zzfc(1,lj1-2))

	endif

!     integrated flux form

	if(fc(lj1-1).ge.0.)then

	  flux(:,lj1-1) = max(0.,xnend(:,1)*zzfc(5,lj1-1)		&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(4,lj1-1)		&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(3,lj1-1)	&
		+ xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(2,lj1-1) 	&
		+ xn_adv(:,(lj1-3)*MAXLIMAX)*zzfc(1,lj1-1))
	  fluxps(lj1-1) = max(0.,psend(1)*zzfc(5,lj1-1)		&
		+ ps3d(lj1*MAXLIMAX)*zzfc(4,lj1-1)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(3,lj1-1)		&
		+ ps3d((lj1-2)*MAXLIMAX)*zzfc(2,lj1-1)		&
		+ ps3d((lj1-3)*MAXLIMAX)*zzfc(1,lj1-1))

	else

	  flux(:,lj1-1) = max(0.,xnend(:,2)*zzfc(5,lj1-1)		&
		+ xnend(:,1)*zzfc(4,lj1-1)			&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(3,lj1-1)		&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(2,lj1-1) 	&
		+ xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(1,lj1-1))
	  fluxps(lj1-1) = max(0.,psend(2)*zzfc(5,lj1-1)		&
		+ psend(1)*zzfc(4,lj1-1)				&
		+ ps3d(lj1*MAXLIMAX)*zzfc(3,lj1-1)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(2,lj1-1)		&
		+ ps3d((lj1-2)*MAXLIMAX)*zzfc(1,lj1-1))

	endif

!     integrated flux form
	if(limthig.eq.lj1)then

	  if(fc(lj1).ge.0.)then

	    flux(:,lj1) = max(0.,xnend(:,2)*zzfc(5,lj1)		&
		+ xnend(:,1)*zzfc(4,lj1)				&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(3,lj1)		&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(2,lj1) 		&
		+ xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(1,lj1))
	    fluxps(lj1) = max(0.,psend(2)*zzfc(5,lj1)		&
		+ psend(1)*zzfc(4,lj1)				&
		+ ps3d(lj1*MAXLIMAX)*zzfc(3,lj1)			&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(2,lj1)		&
		+ ps3d((lj1-2)*MAXLIMAX)*zzfc(1,lj1))

	  else

	    flux(:,lj1) = max(0.,xnend(:,3)*zzfc(5,lj1)	&
		+ xnend(:,2)*zzfc(4,lj1)			&
		+ xnend(:,1)*zzfc(3,lj1)			&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(2,lj1) 	&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1))
	    fluxps(lj1) = max(0.,psend(3)*zzfc(5,lj1)	&
		+ psend(2)*zzfc(4,lj1)			&
		+ psend(1)*zzfc(3,lj1)			&
		+ ps3d(lj1*MAXLIMAX)*zzfc(2,lj1)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1))

	  endif

	else

!     integrated flux form

	  flux(:,lj1) = max(0.,xnend(:,3)*zzfc(5,lj1)	&
		+ xnend(:,2)*zzfc(4,lj1)			&
		+ xnend(:,1)*zzfc(3,lj1)			&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(2,lj1) 	&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1))
	  flux(:,lj1+1) = max(0.,xnend(:,3)*zzfc(5,lj1+1)		&
		+ xnend(:,2)*zzfc(4,lj1+1)			&
		+ xnend(:,1)*zzfc(3,lj1+1)			&
		+ xn_adv(:,lj1*MAXLIMAX)*zzfc(2,lj1+1) 		&
		+ xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1+1))
	  fluxps(lj1) = max(0.,psend(3)*zzfc(5,lj1)	&
		+ psend(2)*zzfc(4,lj1)			&
		+ psend(1)*zzfc(3,lj1)			&
		+ ps3d(lj1*MAXLIMAX)*zzfc(2,lj1)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1))
	  fluxps(lj1+1) = max(0.,psend(3)*zzfc(5,lj1+1)		&
		+ psend(2)*zzfc(4,lj1+1)				&
		+ psend(1)*zzfc(3,lj1+1)				&
		+ ps3d(lj1*MAXLIMAX)*zzfc(2,lj1+1)		&
		+ ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1+1))

	endif

	if(limtlow.eq.-1)then
	  hel1(:) = xnbeg(:,3)*xmdloc(0)
	  hel2(:) = flux(:,0) +  flux(:,-1)
	  where(hel1(:).lt.hel2(:))flux(:,0) 	&
			= flux(:,0)*hel1(:)/hel2(:)
	  hel1ps = psbeg(3)*xmdloc(0)
	  hel2ps = fluxps(0) +  fluxps(-1)
	  if(hel1ps.lt.hel2ps)fluxps(0) 		&
			= fluxps(0)*hel1ps/hel2ps
	  ij = 1
	else
	  if(fc(lj0-1).ge.0.) then
	    flux(:,lj0-1) = amin1(xnbeg(:,3)		&
			*xmdloc(lj0-1)			&
			,flux(:,lj0-1))
	    fluxps(lj0-1) = amin1(psbeg(3)		&
			*xmdloc(lj0-1)			&
			,fluxps(lj0-1))
	    ij = lj0
	  else
	    if(fc(lj0).lt.0.) then
	      flux(:,lj0-1) = -amin1(xn_adv(:,lj0*MAXLIMAX)	&
			*xmdloc(lj0),flux(:,lj0-1))
	      fluxps(lj0-1) = -amin1(ps3d(lj0*MAXLIMAX)		&
			*xmdloc(lj0),fluxps(lj0-1))
	      ij = lj0
	    else
	      hel1(:) = xn_adv(:,lj0*MAXLIMAX)*xmdloc(lj0)
	      hel2(:) = flux(:,lj0) +  flux(:,lj0-1)
	      where(hel1(:).lt.hel2(:))
	        flux(:,lj0-1) = - flux(:,lj0-1)*hel1(:)/hel2(:)
	        flux(:,lj0) = flux(:,lj0)*hel1(:)/hel2(:)
	        xn_adv(:,lj0*MAXLIMAX) = 0.
	      elsewhere
	        flux(:,lj0-1) = -flux(:,lj0-1)
	        xn_adv(:,lj0*MAXLIMAX) =xm2loc(lj0)	&
			*(hel1(:)-hel2(:))
	      end where
	      hel1ps = ps3d(lj0*MAXLIMAX)*xmdloc(lj0)
	      hel2ps = fluxps(lj0) +  fluxps(lj0-1)
	      if(hel1ps.lt.hel2ps)then
	        fluxps(lj0-1) = - fluxps(lj0-1)*hel1ps/hel2ps
	        fluxps(lj0) = fluxps(lj0)*hel1ps/hel2ps
	        ps3d(lj0*MAXLIMAX) = 0.
	      else
	        fluxps(lj0-1) = -fluxps(lj0-1)
	        ps3d(lj0*MAXLIMAX) =xm2loc(lj0)*(hel1ps-hel2ps)
	      endif
	      ij = lj0+1
	    endif
	  endif
	endif

	ijpasses = 0
	do while(.true.)

	  ijpasses = ijpasses+1
	  ijb1(ijpasses) = ij
	  ije1(ijpasses) = -5
	  do while(fc(ij).ge.0.)
	    ije1(ijpasses) = ij
	    ij = ij+1
	    if(ij.gt.lj1-1)then
	      ijb2(ijpasses) = ij
	      ije2(ijpasses) = -5
	      ijb3(ijpasses) = -5
	      goto 257
	    endif
	  enddo
	  ijb2(ijpasses) = ij
	  ije2(ijpasses) = -5
	  do while(fc(ij+1).lt.0.)
	    ije2(ijpasses) = ij
	    ij = ij+1
	    if(ij.gt.lj1-1)then
	      ijb3(ijpasses) = -5
	      goto 257
	    endif
	  enddo
	  ijb3(ijpasses) = ij
	  ij = ij+2
	  if(ij.gt.lj1-1)goto 257
	enddo

257	continue
	ijdoend = .false.
	if(ij.eq.lj1)ijdoend=.true.

	do ijll = 1,ijpasses

	  do ij = ijb1(ijll),ije1(ijll)
	    flux(:,ij) = amin1(xn_adv(:,ij*MAXLIMAX)*xmdloc(ij)	&
			,flux(:,ij))
	    xn_adv(:,ij*MAXLIMAX) =					&
			amax1(0.,xn_adv(:,ij*MAXLIMAX)			&
			-xm2loc(ij)					&
			*(flux(:,ij) - flux(:,ij-1)))
	    fluxps(ij) = amin1(ps3d(ij*MAXLIMAX)*xmdloc(ij)		&
			,fluxps(ij))
	    ps3d(ij*MAXLIMAX) =						&
			amax1(0.,ps3d(ij*MAXLIMAX)			&
			-xm2loc(ij)					&
			*(fluxps(ij) - fluxps(ij-1)))
	  enddo
	  do ij = ijb2(ijll),ije2(ijll)
	    flux(:,ij) = -amin1(xn_adv(:,(ij+1)*MAXLIMAX)*xmdloc(ij+1)&
			,flux(:,ij))
	    xn_adv(:,ij*MAXLIMAX) =					&
			amax1(0.,xn_adv(:,ij*MAXLIMAX)			&
			-xm2loc(ij)*(flux(:,ij)				&
			- flux(:,ij-1)))
	    fluxps(ij) = -amin1(ps3d((ij+1)*MAXLIMAX)*xmdloc(ij+1)	&
			,fluxps(ij))
	    ps3d(ij*MAXLIMAX) =					&
			amax1(0.,ps3d(ij*MAXLIMAX)			&
			-xm2loc(ij)*(fluxps(ij)				&
			- fluxps(ij-1)))
	  enddo
	  ij = ijb3(ijll)
	  if(ij.lt.-3) goto 357
	  hel1(:) = xn_adv(:,(ij+1)*MAXLIMAX)*xmdloc(ij+1)
	  hel2(:) = flux(:,ij+1) +  flux(:,ij)
	  where(hel1(:).lt.hel2(:))
	    flux(:,ij) = - flux(:,ij)*hel1(:)/hel2(:)
	    flux(:,ij+1) = flux(:,ij+1)*hel1(:)/hel2(:)
	    xn_adv(:,(ij+1)*MAXLIMAX) = 0.
	  else where
	    flux(:,ij) = -flux(:,ij)
	    xn_adv(:,(ij+1)*MAXLIMAX) = xm2loc(ij+1)*(hel1(:)-hel2(:))
	  end where
	  xn_adv(:,ij*MAXLIMAX) =				&
			amax1(0.,xn_adv(:,ij*MAXLIMAX)		&
			-xm2loc(ij)*(flux(:,ij) - flux(:,ij-1)))
	  hel1ps = ps3d((ij+1)*MAXLIMAX)*xmdloc(ij+1)
	  hel2ps = fluxps(ij+1) +  fluxps(ij)
	  if(hel1ps.lt.hel2ps)then
	    fluxps(ij) = -fluxps(ij)*hel1ps/hel2ps
	    fluxps(ij+1) = fluxps(ij+1)*hel1ps/hel2ps
	    ps3d((ij+1)*MAXLIMAX) = 0.
	  else
	    fluxps(ij) = -fluxps(ij)
	    ps3d((ij+1)*MAXLIMAX) = xm2loc(ij+1)*(hel1ps-hel2ps)
	  endif
	  ps3d(ij*MAXLIMAX) =					&
			amax1(0.,ps3d(ij*MAXLIMAX)		&
			-xm2loc(ij)*(fluxps(ij) - fluxps(ij-1)))
	enddo

357	continue

	if(ijdoend)then
	  if(limthig.eq.lj1+1)then

	    hel1(:) = xnend(:,1)*xmdloc(lj1+1)
	    hel2(:) = flux(:,lj1+1) + flux(:,lj1)
	    where(hel1(:).lt.hel2(:))
	      flux(:,lj1) 					&
			= -flux(:,lj1)*hel1(:)/hel2(:)
	    else where
	      flux(:,lj1) = -flux(:,lj1)
	    end where
	    xn_adv(:,lj1*MAXLIMAX) =amax1(0.		&
			,xn_adv(:,lj1*MAXLIMAX)-xm2loc(lj1)	&
			*(flux(:,lj1)- flux(:,lj1-1)))
	    hel1ps = psend(1)*xmdloc(lj1+1)
	    hel2ps = fluxps(lj1+1) + fluxps(lj1)
	    if(hel1ps.lt.hel2ps)then
	      fluxps(lj1) 					&
			= -fluxps(lj1)*hel1ps/hel2ps
	    else
	      fluxps(lj1) = -fluxps(lj1)
	    endif
	    ps3d(lj1*MAXLIMAX) =amax1(0.			&
			,ps3d(lj1*MAXLIMAX)-xm2loc(lj1)		&
			*(fluxps(lj1)- fluxps(lj1-1)))

	  else

	    if(fc(lj1).ge.0.) then
	      flux(:,lj1) = amin1(xn_adv(:,lj1*MAXLIMAX)	&
			*xmdloc(lj1),flux(:,lj1))
	      xn_adv(:,lj1*MAXLIMAX) =amax1(0.			&
			,xn_adv(:,lj1*MAXLIMAX)			&
			-xm2loc(lj1)*(flux(:,lj1)- flux(:,lj1-1)))
	      fluxps(lj1) = amin1(ps3d(lj1*MAXLIMAX)		&
			*xmdloc(lj1),fluxps(lj1))
	      ps3d(lj1*MAXLIMAX) =amax1(0.			&
			,ps3d(lj1*MAXLIMAX)			&
			-xm2loc(lj1)*(fluxps(lj1)- fluxps(lj1-1)))
	    else
	      flux(:,lj1) = -amin1(xnend(:,1)*xmdloc(lj1+1)	&
			,flux(:,lj1))
	      xn_adv(:,lj1*MAXLIMAX) =amax1(0.		&
			,xn_adv(:,lj1*MAXLIMAX)			&
			-xm2loc(lj1)				&
			*(flux(:,lj1)- flux(:,lj1-1)))
	      fluxps(lj1) = -amin1(psend(1)*xmdloc(lj1+1)	&
			,fluxps(lj1))
	      ps3d(lj1*MAXLIMAX) =amax1(0.			&
			,ps3d(lj1*MAXLIMAX)-xm2loc(lj1)		&
			*(fluxps(lj1)- fluxps(lj1-1)))
	    endif
	  endif
	endif

!     accumulation of the boundary fluxes

	if (lj0.eq.2) then
	  if(fc(1).ge.0.)then
            fluxin(:) = fluxin(:) + flux(:,1)*fac1
	  else
            fluxout(:) = fluxout(:) - flux(:,1)*fac1
	  endif
	endif

	if (lj1.eq.ljmax-1) then
	  if(fc(lj1).ge.0.)then
            fluxout(:) = fluxout(:) + flux(:,lj1)*fac1
	  else
            fluxin(:) = fluxin(:) - flux(:,lj1)*fac1
	  endif
	endif

	end subroutine advy

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine preadvx(msgnr		&
			,xn_adv,ps3d,vel	&
 			,xnbeg, xnend		&
			,psbeg, psend)

	use Par_ml , only : lj0,lj1,li1		&
			,neighbor,WEST,EAST
	use GenSpec_adv_ml , only : NSPEC_ADV
	implicit none

!	input
	integer,intent(in):: msgnr
	real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX:MAXLIMAX*MAXLJMAX)
	real,intent(in):: ps3d(MAXLIMAX:MAXLIMAX*MAXLJMAX)	&
		,vel(MAXLIMAX+1:(MAXLIMAX+1)*MAXLJMAX)

!	output
	real,intent(out):: xnend(NSPEC_ADV,3,MAXLJMAX)		&
		, xnbeg(NSPEC_ADV,3,MAXLJMAX)
	real,intent(out):: psend(3,MAXLJMAX)			&
		, psbeg(3,MAXLJMAX)

!	local
	integer  i, info

!     Initialize arrays holding boundary slices

!     send to WEST neighbor if any

	if (neighbor(WEST).ge.0) then
	  do i = lj0,lj1
 
	    xnbeg(:,1,i) = xn_adv(:,i*MAXLIMAX)
	    xnbeg(:,2,i) = xn_adv(:,i*MAXLIMAX+1)
	    xnbeg(:,3,i) = xn_adv(:,i*MAXLIMAX+2)

	    psbeg(1,i) = ps3d(i*MAXLIMAX)
	    psbeg(2,i) = ps3d(i*MAXLIMAX+1)
	    psbeg(3,i) = ps3d(i*MAXLIMAX+2)

	  enddo
	  call gc_rsend(msgnr,3*MAXLJMAX*NSPEC_ADV,		&
			neighbor(WEST), info, xnbeg, xnbeg)
	  call gc_rsend(msgnr+100,3*MAXLJMAX,			&
			neighbor(WEST), info, psbeg, psbeg)
	endif

	if (neighbor(EAST).ge.0) then

	  do i = lj0,lj1
	    xnend(:,1,i) = xn_adv(:,i*MAXLIMAX+li1-3)
	    xnend(:,2,i) = xn_adv(:,i*MAXLIMAX+li1-2)
	    xnend(:,3,i) = xn_adv(:,i*MAXLIMAX+li1-1)

	    psend(1,i) = ps3d(i*MAXLIMAX+li1-3)
	    psend(2,i) = ps3d(i*MAXLIMAX+li1-2)
	    psend(3,i) = ps3d(i*MAXLIMAX+li1-1)
	  enddo
	  call gc_rsend(msgnr,3*MAXLJMAX*NSPEC_ADV,		&
			neighbor(EAST), info, xnend, xnend)
	  call gc_rsend(msgnr+100,3*MAXLJMAX,			&
			neighbor(EAST), info, psend, psend)
	endif

	if (neighbor(EAST).lt.0) then

	  do i = lj0,lj1
 	      

	    if(vel(i*(MAXLIMAX+1)+li1).ge.0)then
	      xnend(:,1,i) = 2.*xn_adv(:,i*MAXLIMAX+li1-1)	&
			-xn_adv(:,i*MAXLIMAX+li1-2)
	      xnend(:,2,i) = 3.*xn_adv(:,i*MAXLIMAX+li1-1)	&
			-2.*xn_adv(:,i*MAXLIMAX+li1-2)

	      psend(1,i) = 2.*ps3d(i*MAXLIMAX+li1-1)	&
			-ps3d(i*MAXLIMAX+li1-2)
	      psend(2,i) = 3.*ps3d(i*MAXLIMAX+li1-1)	&
			-2.*ps3d(i*MAXLIMAX+li1-2)
	    else
	      xnend(:,1,i) = xn_adv(:,i*MAXLIMAX+li1)
	      xnend(:,2,i) = xn_adv(:,i*MAXLIMAX+li1)
	      xnend(:,3,i) = xn_adv(:,i*MAXLIMAX+li1)

	      psend(1,i) = ps3d(i*MAXLIMAX+li1)
	      psend(2,i) = ps3d(i*MAXLIMAX+li1)
	      psend(3,i) = ps3d(i*MAXLIMAX+li1)
	    endif
	  enddo
	endif

	if (neighbor(WEST).lt.0) then

	  do i = lj0,lj1
 
	    if(vel(i*(MAXLIMAX+1)+1).lt.0)then
	      xnbeg(:,2,i) = 3.*xn_adv(:,i*MAXLIMAX+1)	&
			-2.*xn_adv(:,i*MAXLIMAX+2)
	      xnbeg(:,3,i) = 2.*xn_adv(:,i*MAXLIMAX+1)	&
			-xn_adv(:,i*MAXLIMAX+2)

	      psbeg(2,i) = 3.*ps3d(i*MAXLIMAX+1)-2.*ps3d(i*MAXLIMAX+2)
	      psbeg(3,i) = 2.*ps3d(i*MAXLIMAX+1)-ps3d(i*MAXLIMAX+2)
	    else

	      xnbeg(:,1,i) = xn_adv(:,i*MAXLIMAX)
	      xnbeg(:,2,i) = xn_adv(:,i*MAXLIMAX)
	      xnbeg(:,3,i) = xn_adv(:,i*MAXLIMAX)
 
	      psbeg(1,i) = ps3d(i*MAXLIMAX)
	      psbeg(2,i) = ps3d(i*MAXLIMAX)
	      psbeg(3,i) = ps3d(i*MAXLIMAX)

	    endif
	  enddo
	endif

	if (neighbor(EAST).ge.0) then
	  call gc_rrecv(msgnr,MAXLJMAX*3*NSPEC_ADV,		&
			neighbor(EAST), info, xnend, xnend)
	  call gc_rrecv(msgnr+100,MAXLJMAX*3,			&
			neighbor(EAST), info, psend, psend)
	endif

!     receive from WEST neighbor if any

	if (neighbor(WEST).ge.0) then
	  call gc_rrecv(msgnr,MAXLJMAX*3*NSPEC_ADV,		&
			neighbor(WEST), info, xnbeg, xnbeg)
	  call gc_rrecv(msgnr+100,MAXLJMAX*3,			&
			neighbor(WEST), info, psbeg, psbeg)
	endif

	end subroutine preadvx

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine preadvy(msgnr		&
			,xn_adv,ps3d,vel	&
			,xnbeg, xnend		&
			,psbeg, psend)

	use Par_ml , only : li0,li1,lj1		&
			,neighbor,NORTH,SOUTH
	use GenSpec_adv_ml , only : NSPEC_ADV
	implicit none

!	input
	integer,intent(in):: msgnr
	real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX*MAXLJMAX)
	real,intent(in):: ps3d(MAXLIMAX*MAXLJMAX)		&
		,vel(MAXLIMAX*(MAXLJMAX+1))

!	output
	real,intent(out):: xnend(NSPEC_ADV,3,MAXLIMAX)		&
		, xnbeg(NSPEC_ADV,3,MAXLIMAX)
	real,intent(out):: psend(3,MAXLIMAX)			&		
		, psbeg(3,MAXLIMAX)

!	local
	integer  i, info

!     Initialize arrays holding boundary slices

!     send to SOUTH neighbor if any

	if (neighbor(SOUTH).ge.0) then
	  do i = li0,li1
 
	    xnbeg(:,1,i) = xn_adv(:,i)
	    xnbeg(:,2,i) = xn_adv(:,i+MAXLIMAX)
	    xnbeg(:,3,i) = xn_adv(:,i+2*MAXLIMAX)
 
	    psbeg(1,i) = ps3d(i)
	    psbeg(2,i) = ps3d(i+MAXLIMAX)
	    psbeg(3,i) = ps3d(i+2*MAXLIMAX)

	  enddo
	  call gc_rsend(msgnr,3*MAXLIMAX*NSPEC_ADV,		&
			neighbor(SOUTH), info, xnbeg, xnbeg)
	  call gc_rsend(msgnr+100,3*MAXLIMAX,			&
			neighbor(SOUTH), info, psbeg, psbeg)

	endif

	if (neighbor(NORTH).ge.0) then

	  do i = li0,li1
	    xnend(:,1,i) = xn_adv(:,i+(lj1-3)*MAXLIMAX)
	    xnend(:,2,i) = xn_adv(:,i+(lj1-2)*MAXLIMAX)
	    xnend(:,3,i) = xn_adv(:,i+(lj1-1)*MAXLIMAX)

	    psend(1,i) = ps3d(i+(lj1-3)*MAXLIMAX)
	    psend(2,i) = ps3d(i+(lj1-2)*MAXLIMAX)
	    psend(3,i) = ps3d(i+(lj1-1)*MAXLIMAX)
	  enddo
	  call gc_rsend(msgnr,3*MAXLIMAX*NSPEC_ADV,		&  
			neighbor(NORTH), info, xnend, xnend)
	  call gc_rsend(msgnr+100,3*MAXLIMAX,			&
			neighbor(NORTH), info, psend, psend)
	endif

	if (neighbor(SOUTH).lt.0) then

	  do i = li0,li1

	    if(vel(i+MAXLIMAX).lt.0)then
	      xnbeg(:,2,i) = 3.*xn_adv(:,i+MAXLIMAX)	&
			-2.*xn_adv(:,i+2*MAXLIMAX)
	      xnbeg(:,3,i) = 2.*xn_adv(:,i+MAXLIMAX)	&
			-xn_adv(:,i+2*MAXLIMAX)

	      psbeg(2,i) = 3.*ps3d(i+MAXLIMAX)-2.*ps3d(i+2*MAXLIMAX)
	      psbeg(3,i) = 2.*ps3d(i+MAXLIMAX)-ps3d(i+2*MAXLIMAX)
	    else
 
	      xnbeg(:,1,i) = xn_adv(:,i)
	      xnbeg(:,2,i) = xn_adv(:,i)
	      xnbeg(:,3,i) = xn_adv(:,i)
 
	      psbeg(1,i) = ps3d(i)
	      psbeg(2,i) = ps3d(i)
	      psbeg(3,i) = ps3d(i)

	    endif
	  enddo
	endif

	if (neighbor(NORTH).lt.0) then

	  do i = li0,li1
 
	    if(vel(i+lj1*MAXLIMAX).ge.0)then
	      xnend(:,1,i) = 2.*xn_adv(:,i+(lj1-1)*MAXLIMAX)	&
			-xn_adv(:,i+(lj1-2)*MAXLIMAX)
	      xnend(:,2,i) = 3.*xn_adv(:,i+(lj1-1)*MAXLIMAX)	&
			-2.*xn_adv(:,i+(lj1-2)*MAXLIMAX)

	      psend(1,i) = 2.*ps3d(i+(lj1-1)*MAXLIMAX)	&
			-ps3d(i+(lj1-2)*MAXLIMAX)
	      psend(2,i) = 3.*ps3d(i+(lj1-1)*MAXLIMAX)	&
			-2.*ps3d(i+(lj1-2)*MAXLIMAX)
	    else
	      xnend(:,1,i) = xn_adv(:,i+lj1*MAXLIMAX)
	      xnend(:,2,i) = xn_adv(:,i+lj1*MAXLIMAX)
	      xnend(:,3,i) = xn_adv(:,i+lj1*MAXLIMAX)

	      psend(1,i) = ps3d(i+lj1*MAXLIMAX)
	      psend(2,i) = ps3d(i+lj1*MAXLIMAX)
	      psend(3,i) = ps3d(i+lj1*MAXLIMAX)
	    endif
	  enddo
	endif

	if (neighbor(NORTH).ge.0) then
	  call gc_rrecv(msgnr,MAXLIMAX*3*NSPEC_ADV,		&
			neighbor(NORTH), info, xnend, xnend)
	  call gc_rrecv(msgnr+100,MAXLIMAX*3,			&
			neighbor(NORTH), info, psend, psend)

	endif

!     receive from SOUTH neighbor if any

	if (neighbor(SOUTH).ge.0) then
	  call gc_rrecv(msgnr,MAXLIMAX*3*NSPEC_ADV,		&
			neighbor(SOUTH), info, xnbeg, xnbeg)
	  call gc_rrecv(msgnr+100,MAXLIMAX*3,			&
			neighbor(SOUTH), info, psbeg, psbeg)
	endif

	end subroutine preadvy

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine adv_var(numt)
!hf u2 	use My_Runmode_ml, only : stop_test
	use Par_ml , only : limax,ljmax,NPROC,me	&
		,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC			&
		,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
	use GridValues_ml , only : GRIDWIDTH_M
	use Met_ml  , only : u,v,sdot
	use ModelConstants_ml, only : dt_advec

	integer,intent(in)::  numt

!	local

	integer i,j,k,info,nr
	real dhskmax,sdotmax,sdotmin
	real sdotmaxk,sdotmink
	real ulmin,ulmax,vlmin,vlmax

	nr = 2
	if (numt.eq.1) nr = 1

!     send to WEST neighbor if any

	if (neighbor(WEST) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do j = 1,ljmax
	      uw(j,k,nr) = u(1,j,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_EAST2, MAXLJMAX*KMAX_MID, neighbor(WEST),      &
		info, ue(1,1,nr), uw(1,1,nr))
	endif

!     send to EAST neighbor if any

	if (neighbor(EAST) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do j = 1,ljmax
	      ue(j,k,nr) = u(limax-1,j,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_WEST2, MAXLJMAX*KMAX_MID, neighbor(EAST),      &
		info, uw(1,1,nr), ue(1,1,nr))
	endif

!     send to SOUTH neighbor if any

	if (neighbor(SOUTH) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do i = 1,limax
	      vs(i,k,nr) = v(i,1,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_NORTH2, MAXLIMAX*KMAX_MID, neighbor(SOUTH),      &
		info, vn(1,1,nr), vs(1,1,nr))
	endif

!     send to NORTH neighbor if any

	if (neighbor(NORTH) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do i = 1,limax
	      vn(i,k,nr) = v(i,ljmax-1,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_SOUTH2, MAXLIMAX*KMAX_MID, neighbor(NORTH),      &
		info, vs(1,1,nr), vn(1,1,nr))
	endif

!     receive from EAST neighbor if any

	if (neighbor(EAST) .ne. NOPROC) then
        call gc_rrecv(MSG_EAST2, MAXLJMAX*KMAX_MID, neighbor(EAST),      &
		info, ue(1,1,nr), uw(1,1,nr))
	endif

!     receive from WEST neighbor if any

	if (neighbor(WEST) .ne. NOPROC) then
        call gc_rrecv(MSG_WEST2, MAXLJMAX*KMAX_MID, neighbor(WEST),      &
		info, uw(1,1,nr), ue(1,1,nr))
	endif

!     receive from NORTH neighbor if any

	if (neighbor(NORTH) .ne. NOPROC) then
        call gc_rrecv(MSG_NORTH2, MAXLIMAX*KMAX_MID, neighbor(NORTH),      &
		info, vn(1,1,nr), vs(1,1,nr))
	endif

!     receive from SOUTH neighbor if any

	if (neighbor(SOUTH) .ne. NOPROC) then
        call gc_rrecv(MSG_SOUTH2, MAXLIMAX*KMAX_MID, neighbor(SOUTH),      &
		info, vs(1,1,nr), vn(1,1,nr))
	endif

	return


	end subroutine adv_var
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine adv_int

	use ModelConstants_ml , only : nmax,nstep
	implicit none

	real div

	if (nstep.lt.nmax) then
	  div = 1./real(nmax-(nstep-1))
	  ue(:,:,1) = ue(:,:,1) + (ue(:,:,2) - ue(:,:,1))*div
	  uw(:,:,1) = uw(:,:,1) + (uw(:,:,2) - uw(:,:,1))*div
	  vs(:,:,1) = vs(:,:,1) + (vs(:,:,2) - vs(:,:,1))*div
	  vn(:,:,1) = vn(:,:,1) + (vn(:,:,2) - vn(:,:,1))*div

	else

	  ue(:,:,1) = ue(:,:,2)
	  uw(:,:,1) = uw(:,:,2)
	  vs(:,:,1) = vs(:,:,2)
	  vn(:,:,1) = vn(:,:,2)

	endif

	end subroutine adv_int
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Advection_ml
