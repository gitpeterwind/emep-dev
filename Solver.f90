module Chemsolver_ml
!
! modified by Peter february 2003:
! variable timestep
! decoupling of (NO3,N2O5), (PAN,CH3COO2), (MPAN,MACRO2)
! removed xold=... when(xnew>CPINIT)
!


  implicit none
  private

  public  :: chemistry        ! Runs chemical solver

integer, parameter::nchemMAX=12

contains
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
  subroutine chemistry(i,j,dt,Niter,reset_chem)
    !hf MADE
    use Radiation_ml,          only :  zen ! zenith angle, degrees
    use GenSpec_tot_ml     ! => NSPEC_TOT, O3, NO2, etc.
    !u2 !hf MADE
    use GenSpec_bgn_ml      ! => IXBGN_  indices and xn_2d_bgn values
    use GenRates_rct_ml,    only : set_night_rct, ONLY_NIGHT
    use Par_ml,              only : me,MAXLIMAX,MAXLJMAX,li1,lj1,li0,lj0  ! me for TEST
    use ModelConstants_ml,   only : KMAX_MID, KCHEMTOP ,dt_advec
    use Setup_1dfields_ml,   only : &
         rcemis,izen       & ! photolysis, emissions
         ,rcbio                & ! biogenic emis
         !u1 ,rct, rctroe, rcmisc, & ! rate-coeffients
    ,rct, rcmisc, & ! rate-coeffients
         xn_2d,&               ! rename xn_2d to use simply x inside chemistry
         !hf u2
    rh
    !u2 !hf MADE 
    !u2          xn_2d_bgn
    use Aqueous_ml,        only : &
         aqrck,                    &
         ICLOHSO2  &
         ,ICLRC1   ,ICLRC2   ,ICLRC3   
    use Biogenics_ml, only: BIO_ISOP, BIO_TERP

    use My_Emis_ml                        ! => QRCNO, etc.
    use DefPhotolysis_ml                  ! => IDHNO3, etc.
    use Emissions_ml, only : KEMISTOP     !rv1.2.1 change
    use OrganicAerosol_ml, only : Fgas
    !
    !**********************************************************************
    ! The following solver is a simplified version of that suggested by
    ! Jan Verwer. First, we assume fixed time steps,  and following
    ! Frode Flatoy we use Dchem to keep track of changes from call to call.
    ! Steffen has pointed out that fewer variables are needed during
    ! the actual integration. In order to compare with previous EMEP versions
    ! and Jan Verwer's equations, here are the two methods written side-
    ! by-side:
    !
    !========= su ==================== common ============== orig-emep ========!
    !                             xold = x - Dchem_2d   
    !                             xold = max( xold ,0.0)    
    !                     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                     C    Start of integration loop      C
    !                     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !                Calculates the initial iterative by extrapolation:
    ! xold = ( 4.0*x - xold) * inv3         |
    ! xnew = 3.0 * xold - 2.0 * x           |  xnew = 2 * x - xold
    !                            where ( xnew < CPINIT  ) 
    !                               xnew = CPINIT
    ! xold = ( 2 * x - CPINIT ) * inv3      !  xold = 2 * x - CPINIT   
    !
    !                         GenChem output equations here:
    !                           do iter = 1, Numiter
    !                       include 'My_Reactions.inc'
    !                                EQUATIONS ARE
    ! xnew = xold+dt2*P                     |   xnew = (4*x-xold)*2.dt.P
    !        ----------                     |          -----------------
    !        1+dt2*L                        |            3+2dt.L
    !                           end do ! End iterations
    !
    !                            if ( ichem <  nchem ) then   
    ! xold = ( 4.0 * xnew - x) * inv3       !  xold = x
    !                                 x    = max( xnew, 0.0 )
    !                            end if
    !                       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                       C    End of integration loop      C
    !                       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                 Saves Dchem and returns the new concentrations:
    !                             Dchem_2d = xnew   - x
    !fix- may not be needed, but... x        = xnew
    !
    !========= su ========================================== jw ===============!
    !
    !.. In
    integer, intent(in) ::  i,j       ! Coordinates (needed for Dchem)
    integer, intent(in) ::  dt        ! Local timestep in chemistry (=dt_chem)
    integer, intent(in) ::  Niter     !  No. iterations used in 2step
    logical, intent(in) ::  reset_chem ! Increases iterations, e.g. at start

    !su
    real, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), save :: Dchem=0.

    logical, save ::  first_call = .true.

    real, parameter ::  CPINIT = 0.0 ! 1.0e-30  ! small value for init
    real, parameter :: dt_advec_inv=1./dt_advec

    !  Local
    integer ::  k, ns, ichem, iter,n   ! Loop indices
    !  integer ::  Numiter              ! Loop index for iterations
    integer,save ::  nchem                ! No chem time-steps
    real    ::  dt2    ,dt2save              ! For su2step, 2/3* dt
    real    ::  inv3                  ! For su2step, 1.0/3.0
    real    ::  P, L                ! Production, loss terms
    ! su's suggestion to help include My_Reactions only once:
    integer, dimension(KCHEMTOP:KMAX_MID) :: toiter

    real LVOC
    real LNO,LSO2,LC2H5,LHCHO,LMARO,LMAL
    real PASO,PSOA1,PSOA3,PCH4
    real hct13no,hccho2,hct3oh
    real hct30c2h4,hct31c3h6,hct32iso,hct37mvk
    real hpsoa1,hpsoa3,hccoxyl,hctisomvk
    real hccmglyox,hphomglyox,hccglyox,hphoglyox
    real hphoch3cho,hct24ch3cho,hct12ch4,hct28c2h5oh
    real hrct28oh,hrct16ch3oh,hrct16oh,hrct14ch3o2,hrct15ch3o2
    real POXYO2
    real xextrapol,L1,L2,P1,P2,C1,C2,DIVID !help variable
!    real Ntotnew,Ntotold,Nrel,NOy,NOynew
!    real,save :: Nrelmin=1.E40,Nrelmax=0.

    ! Concentrations : xold=old, x=current, xnew=predicted
    ! - dimensioned to have same size as "x" 

    real, dimension(NSPEC_TOT) ::  x, xold ,xnew     ! Working array [molecules/cm3]
    real, dimension(nchemMAX), save :: dti ! variable timestep*(c+1)/(c+2)
    real, dimension(nchemMAX), save :: coeff1,coeff2,cc ! coefficients for variable timestep

    integer :: nextraiter

    if ( first_call ) then
       call makedt(dti,nchem,coeff1,coeff2,cc,dt_advec)
       first_call = .false.
    endif
    !su2step
    inv3 = 1.0/3.0
    dt2  = 2.0 * inv3 * dt

    ! su's suggetsion to help include My_Reactions only once:
    !  toiter gives the number of iterations used in TWOSTEP. Use
    !  more iterations near ground:
    !rv1.2.1 - changed to KEMISTOP

    toiter(KCHEMTOP:5)    = 1    ! Upper levels - slow chemistry
    toiter(6:KEMISTOP-1)     = 2    ! Medium and cloud levels 
    toiter(KEMISTOP:KMAX_MID) = 3    !  Near-ground, emis levels

    !...  if we run with dt_chem less than dt_advec sec., we start a loop here.

    !pw    nchem = nint(dt_advec)/dt

    !       Numiter = Niter 
    !       if( reset_chem ) then
    !           Numiter = max(3,Niter) 
    !       end if

    !hf u2 
    !**/ Only NO2+O3->H+ +NO3- at night time and in the
    !    8 lowest layers and if rh>0.5
    if (ONLY_NIGHT) call set_night_rct(rct,rh,i,j)


    ! Establishment of initial conditions:
    ! Previous concentrations are estimated by the current
    ! minus Dchem because the current may be changed by
    ! processes outside the chemistry:

    do k = 2, KMAX_MID
       ! Remember, xnew   (:,:) =~ xn_2d(:,:)

       xnew(:) = xn_2d(:,k)

       !      x(:) = xn_2d(:,k) - Dchem(:,k,i,j)
       x(:) = xn_2d(:,k) - Dchem(:,k,i,j)*dti(1)*1.5
       x(:) = max (x(:), 0.0)

       !        if ( first_call .and. me == 0 ) then
       !            print *, "TTTT SOLVER HEAD  ", first_call
       !            print *, "TTTT SOLVER CHEMSIZE  ", CHEMSIZE
       !            print *, "TTTT SOLVER HEAD H2, x  ", H2,  x(H2 ,10)
       !            print *, "TTTT SOLVER HEAD H2 old ", H2,  xold(H2 ,10)
       !            print *, "TTTT SOLVER HEAD H2 D_  ", H2,  Dchem_2d(H2 ,10)
       !            print *, "TTTT SOLVER HEAD O3, x  ", O3,  x(O3 ,10)
       !            print *, "TTTT SOLVER HEAD O3 old ", O3,  xold(O3 ,10)
       !            print *, "TTTT SOLVER HEAD O3 D_  ", O3,  Dchem_2d(O3 ,10)
       !            print *, "TTTT size 1,2 = ", size(x,1), size(x,2)
       !        end if



       !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       !     Start of integration loop      C
       !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       do ichem = 1, nchem

          !         xold(:) = ( 4.0*xnew(:) - x(:) ) * inv3
          !         x(:)    = xnew(:)
          !       xnew(:) = 3.0 * xold(:) - 2.0 * x(:)
          do n=1,NSPEC_TOT
             xextrapol = xnew(n) + (xnew(n)-x(n)) *cc(ichem)
             xold(n) = coeff1(ichem)*xnew(n) - coeff2(ichem)*x(n) 
             x(n) = xnew(n)
             xnew(n) = xextrapol
          enddo

          !      dt2  = 2.0 * inv3 * dt
          dt2  =  dti(ichem) !*(1.+cc(ichem))/(1.+2.*cc(ichem))

          where ( xnew(:) < CPINIT  ) 
             xnew(:) = CPINIT
!             !su2step xold(:,:) = 2 * x(:,:) - CPINIT   
!             !           xold = ( 2 * x - CPINIT ) * inv3   ????
!             !           xold(:) = ( 2 * x(:) + CPINIT ) * inv3   !pw why change the past?
!             !           xold(:) = coeff1(ichem)*x(:)-coeff2(ichem)*( x(:) + (x(:)-xnew(:))/cc(ichem) )!pw why change the past?
          end where

          include 'My_Reactions.inc' 

       end do ! Main loop
       !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       !     End of integration loop      C
       !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       !  Saves Dchem and returns the new concentrations:

       !pw:use rate of change instead of asbolute changes   Dchem(:,k,i,j) = xnew(:)   - x(:)

       Dchem(:,k,i,j) = (xnew(:)   - xn_2d(:,k))*dt_advec_inv
       xn_2d(:,k) = xnew(:)

       !>>>>>>>>>>>>>>>>>>>>>>>>>>>

    enddo !k

    !        if ( first_call .and. me == 0 ) then
    !            print *, "TTTT END OF SOLVER ", first_call
    !            print *, "TTTT SOLVER H2     ", H2,  x   (H2 ,10)
    !            print *, "TTTT SOLVER H2 old ", H2,  xold(H2 ,10)
    !            print *, "TTTT SOLVER H2 new ", H2,  xnew(H2 ,10)

    !            print *, "TTTT SOLVER CH4    ", CH4, x(CH4,10)
    !            print *, "TTTT SOLVER O3     ", O3,  x(O3 ,10)
    !          !  first_call = .false.
    !        endif
    !        first_call = .false.
  end subroutine chemistry

 !--------------------------------------------------------------

! function Cnew( dt, P, L, Cn_1, Cn ) result (C_new)
! !--------------------------------------------------------------
! !+   TWOSTEP Interation formula based on a Gauss-Seidel method
! !--------------------------------------------------------------
!      integer, intent(in) ::  dt
!      real, dimension( : ),  intent(in) ::  P, L, Cn_1, Cn
!      real, dimension( size(P) ) ::  C_new
!      C_new = ( 4.0*Cn - Cn_1 + 2.0*dt*P )/ &
!                ( 3.0 + 2.0*dt*L )
!!      if (L.gt.(10./dt)) Cnew = P/L
!      end function Cnew

subroutine  makedt(dti,nchem,coeff1,coeff2,cc,dt_tot)
!
! make coefficients for two-step
! written by Peter February 2003
!
! The formulas for coeff1, coeff2 and dti can be found in:
!  J.G. Verwer and D. Simpson, "Explicit methods for stiff ODEs from
!  atmospheric chemistry", Aplied Numerical Mathematics 18 (1995) 413
!
!Note: It is better to take first some small steps, and then
!larger steps, than to increase the timestep gradually.
!Is this a sign that there is a bug in the coefficients?
!

    use Par_ml,              only : me,NPROC
    use ModelConstants_ml,   only : dt_advec

implicit none

real, dimension(nchemMAX),intent(out) :: dti,coeff1,coeff2,cc
real, intent(in) ::dt_tot
integer, intent(out) ::nchem
real :: ttot,dt_first,dt_max,dtleft,tleft,step,dt(nchemMAX)
integer ::i,j

!dt_first=8.
!dt_max=360.
!nchem=14
!step=exp(log(dt_max/dt_first)/(nchem-1))
!write(*,*)'step ',step
!step=1.4
!nchem=max(1,int(log(dt_max/dt_first)/log(step)))+1
!nchem=12
!

!ttot=0.
!dt(1)=dt_first
!ttot=dt(1)
!tleft=dt_tot-ttot
!dtleft=tleft/(nchem-1)
!do i=2,nchem-1
!   if(dtleft<dt(i-1)*step)then
!      dt(i)=dtleft
!   else
!      dt(i)=dt(i-1)*step
!   endif
!   ttot=ttot+dt(i)
!   tleft=dt_tot-ttot
!   dtleft=tleft/(nchem-i)
!   if(dtleft>dt_max)then
!      ttot=ttot-dt(i)
!      tleft=dt_tot-ttot
!      dtleft=tleft/float(nchem-i+1)
!      dt(i)=dtleft
!      ttot=ttot+dt(i)
!      tleft=dt_tot-ttot
!      dtleft=tleft/(nchem-i)
!   endif  
!enddo
!
!dt(nchem)=dt_tot-ttot
!ttot=ttot+dt(nchem)

nchem=12

!dt=1100./7.
dt=(dt_advec-100.)/(nchem-5)
dt(1)=20.
dt(2)=20.
dt(3)=20.
dt(4)=20.
dt(5)=20.

if(dt_advec<520.)then
   nchem=5+int((dt_advec-100.)/60.)
   dt=(dt_advec-100.)/(nchem-5)
   dt(1)=20.
   dt(2)=20.
   dt(3)=20.
   dt(4)=20.
   dt(5)=20.
endif
if(dt_advec<=100.)then
   nchem=int(dt_advec/20.)+1
   dt=(dt_advec)/(nchem)
endif
if(dt_advec<20.)then
call gc_abort(me,NPROC, "makedt: dt_advec too small ")
endif


if(nchem > nchemMAX)print *,'WARNING: nchemMAX too small'
nchem=min(nchemMAX,nchem)
if(me == 0)then
write(*,*)'Number of timesteps in Solver: ',nchem
27 format('timestep ',I,F13.6,' total: ',F13.6)
ttot=0.0
do i=1,nchem
   ttot=ttot+dt(i)
   write(*,27)i,dt(i),ttot
enddo

!check that we are using consistent timesteps
if(abs(ttot-dt_advec)>1.E-5)then
print *,'dt_advec and dt not compatible'
call gc_abort(me,NPROC, "makedt: error ")
endif
endif

!dt(0)=dt(1)
cc(1)=1.
coeff2(1)=1./(cc(1)**2+2*cc(1))
coeff1(1)=(cc(1)+1)**2*coeff2(1)
dti(1)=((cc(1)+1)/(cc(1)+2))*dt(1)
do i=2,nchem
cc(i)=dt(i-1)/dt(i)
coeff2(i)=1./(cc(i)**2+2*cc(i))
coeff1(i)=((cc(i)+1)**2)*coeff2(i)
dti(i)=((cc(i)+1)/(cc(i)+2))*dt(i)
!dti(i)=dt(i)
cc(i)=1./cc(i)
enddo

end subroutine makedt

end module Chemsolver_ml
