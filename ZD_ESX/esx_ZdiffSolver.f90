! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!> MODULE
!! Solves vertical diffusion in ESX system
!@author
!> Juha-Pekka Tuovinen and David Simpson
!> 2013
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

module esx_ZdiffSolver
  use esx_MassBudget, only: mass, GRNDDEP,  TOPLOSS, CANDEP, TOT, print_mass
  use esx_Variables, only: esx, Zmet
  implicit none
  private

  public  :: ZdiffSolver
  private :: inv_3diag
  private :: get_mass


  ! weight of (j+1)th and jth time step (Crank-Nicolson, if wj1=0.5)
   real, private, parameter ::  wj1 = 0.5, wj = 1-wj1  

   integer, public, save  :: nDiffSteps = 0 ! external use
   integer, private, save :: nSteps     = 0 ! internal use

contains

  subroutine ZdiffSolver(nz,cindex,dti,Vd,Ve,D,E,Fb,Ft,fixedBC,nSubSteps,concn,debug_level)
   integer, intent(in) :: nz      !! Number of layers
   integer, intent(in) :: cindex  !! index of chemical spec, for get_mass
   real, intent(in) :: &
      dti     &! time step wanted [s]
     ,Vd      &! deposition velocity [m/s]
     ,Ve      &! escape velocity [m/s]
   ! Fluxes, should have units that match chem-solver, but per m2. Could be
   ! molecules/m2/s or (for tests) unitless (#/m2/s), hence:
     ,Fb      &! bottom boundary flux [X/(m2s)] (+ = down)
     ,Ft       ! top    boundary flux [X/(m2s)] (+ = up)

   real, dimension(nz), intent(in) :: &
      D        &! 1st order source/sink coefficient [1/s]     = S1
     ,E         ! 0th order source/sink coefficient [g/(m3s)] = S0

   real, dimension(nz), intent(inout) :: concn
   real, dimension(nz) :: old
   logical, intent(in) :: fixedBC  !  fixed top boundary conc
   integer, intent(in) ::  debug_level  ! 0, 1 or 2

   real, dimension(nz) :: &
      R               &! 
     ,Ua, La, Ma      &! 
     ,U               &!upper diagonal,  U(1:nz+1), U(nz+1)=0
     ,L               &!lower diagonal,  L(1:nz+1), L(1)=0
     ,M               !main diagonal,   M(1:nz+1)

   real ::   Vda, Vea, Fba, Fta, stabtest, dt, dtstab
   integer :: j, k, k2,  nz1, nz2
   integer :: nSubSteps  ! number of time steps in C-N method, to be calculated
   logical, save :: first_call = .true.

   real, pointer, dimension(:) :: z, dz, dzmid, Kz
!
   z     => esx%z(1:nz)
   dz    => esx%dz(1:nz)
   dzmid => esx%dzmid(1:nz)
   Kz    => Zmet(1:nz)%Kz

   nz1 = nz -1  ! Number of mid layers, e.g. for  Kz
   nz2 = nz -2  ! Used in rhs calculations

  !> Stability check for ZdiffSolver time-step
  !! To prevent oscillations, we need  K.dt/(dz**2) < 0.5, or dt < 0.5 dz**2/K.
  !
   if ( dt > tiny(dt) ) dt = min( dt, dti )  !> Selects shortest dt
   dtstab = 0.5*minval( dz(1:nz1)**2 / Kz(1:nz1) )
   dt = min( esx%dt_Zdiff, dtstab ) 
   nSubSteps = ceiling( dti/dt )
   nSubSteps = max(1, nSubSteps)
   dt = dti/nSubSteps

   if( first_call .and. debug_level >0 ) then
     print "(a,3g12.3,a,g12.3,i6)", "ZdiffDts:", &
       dti,dtstab,esx%dt_Zdiff, " => dt, nSubSteps", dt, nSubSteps

     if( debug_level>1 ) then !! Show where the stability criteria max is
       print "(a,i6,20es10.3)", "debug-Zdiff start:",nSubSteps, dt
       print "(a,a6,20a10)", "Zdiff start:","k", "%z(k)",  "dzmid(k)",&
          "dz(k)", "D(k)", "E(k)", "K(k)", "concn(k)", "Crit."

       do k =  nz1, 1, -1
         stabtest =  Kz(k) * dti/dz(k)**2
         print "(a,i6,3f8.2,20es10.2)", "Kdiff start:",k, z(k+1), dzmid(k), &
          dz(k), D(k), E(k), Kz(k), concn(k),  stabtest
       end do
     end if
   end if

 !>  start solution

  Ua(1:nz1)   = Kz(1:nz1)*dt/(dzmid(1:nz1)*dz(1:nz1))
  Ua(nz)=0
  La(2:nz)     = Kz(1:nz1)*dt/(dzmid(1:nz1)*dz(2:nz))
  La(1)=0
  Ma(1:nz)     = D(1:nz)*dt
  Vda = Vd*dt/dz(1)
  Vea = Ve*dt/dz(nz)
  Fba = Fb/ dz(1)
  Fta = Ft/ dz(nz)

 !> 3-diagonal matrix
    U = -wj1*Ua             !upper diagonal,  U(1:nz), U(nz)=0
    L = -wj1*La             !lower diagonal,  L(1:nz), L(1)=0
    M = 1 - U - L - wj1*Ma  !main diagonal,   M(1:nz)
    M(1)  = M(1) + wj1*Vda 
    M(nz) = M(nz) + wj1*Vea

    if ( fixedBC ) then
        L(nz) = 0.0 ! SCOPE CHANGED. kdiff had n+1 for both
        M(nz) = 1.0
    end if

 !> ---start time integration---
    do  j = 1, nSubSteps

       !>>>right-hand side:    
       R(1)    = wj*Ua(1)                        * concn(2)      +&
                 (1 - wj*(Ua(1) + Vda - Ma(1)))  * concn(1)      +&
                 ( E(1) -Fba )  *dt
       R(2:nz1) = wj*Ua(2:nz1)                                * concn(3:nz) +&
                 (1 - wj*(Ua(2:nz1) - Ma(2:nz1) + La(2:nz1))) * concn(2:nz1) +&
                 wj*La(2:nz1)                                 * concn(1:nz2) +&
                 E(2:nz1)*dt
       R(nz) = (1 - wj*(La(nz) + Vea - Ma(nz)))  * concn(nz)   +&
                 wj*La(nz)                       * concn(nz1)  +&
                 ( E(nz) - Fta ) *dt
       !<<

       if ( fixedBC ) R(nz) = concn(nz) ! FIXED top boundary concentration

       old(:) = concn(:)

       concn = inv_3diag(U,M,L,R)  !solve c(j+1)

       if( debug_level>1 ) then
        print "(a,f12.3)", "Kdiff-changes ------------------------------ t=" , esx%Time
        do k =1, nz
          print "(a,2i3,20es12.3)", "Kdiff-changes ", j,k,  &
                old(k), concn(k), R(k), Ua(k), E(k), Ma(k)
        end do
       end if

     if ( any( concn < 0.0) ) then
       print *, "NEG!!! ", j
       do k2 =  nz, 1, -1
         k=min(k2,nz-1)
         print "(a,i6,20es10.2)", "Kdiff fail!:",k, z(k+1), dzmid(k), &
        dz(k), D(k), E(k), Kz(k), concn(k),  stabtest
       end do
       stop "NEG!!!"
     end if

     !>>testing:
       ! mass(j) = sum(concn*dz)   !   int_mass(concn,dz);
       !anmass(j) = an_sol_m1(j*dt,z,Ka,Vd);
       !anmass(j) = exp(D(1)*j*dt); %D(z)=D(1), c(z)=const%
       !<<<

       ! Get mass-budget terms for this species and time-step

        nSteps = nDiffSteps+j
        call get_mass(cindex,nz,Vd,Ve,Fb,Ft,old,concn,D,dt,debug_level)

       first_call = .false.

    end do
    !---end time integration---


end subroutine ZdiffSolver
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inv_3diag(U,M,L,R) result(x)
  !x = inv_3diag(U,M,L,R)
  ! U(1:n) = upper diagonal, U(n)=0
  ! M(1:n) = main diagonal
  ! L(1:n) = lower diagonal, L(1)=0
  ! R(1:n) = right-hand side
 
  real, intent(in), dimension(:) :: U, M, L, R
  real, dimension(size(U)) :: x
  real, dimension(size(U)) :: A, B
  real :: div
  integer :: n, i

  n = size(U);
 
  A(1) = -U(1)/M(1)
  B(1) =  R(1)/M(1)
  do  i = 2, n
    div = M(i) + L(i)*A(i-1)
    A(i) = -U(i)/div   !A(n) not used
    B(i) = (R(i) - L(i)*B(i-1))/div
  end do
 
  x(n) = B(n);
  do i = n-1,1, -1 ! QUERY
     x(i) = A(i)*x(i+1) + B(i)
  end do
end function inv_3diag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_mass(ind,nz,Vd,Ve,Fb,Ft,c,c1,D,dt,debug_level)
!!![massVd massVe massD cumD1] = out_mass(w,w1,Vd,Ve,Fb,Ft,c,c1,D,dz,dt,cumD)
integer, intent(in) :: ind, nz     ! species number,levels
integer, intent(in) :: debug_level
real, intent(in) :: Vd,Ve,Fb,Ft,dt
real, intent(in), dimension(:) :: c,c1,D
real, pointer, dimension(:) :: dz

 dz    => esx%dz(1:nz)

 associate ( m=> mass(ind) )

   m%b(TOT)     = sum(c*dz)

   ! Time-step loss terms:

   m%zloss(1:nz) = -D(:) *( wj*c(:) + wj1*c1(:)) *dz(:)*dt! canopy deposition

   m%loss(GRNDDEP) = (Vd*(wj*c(1)  + wj1*c1(1) ) + Fb)*dt 
   m%loss(TOPLOSS) = (Ve*(wj*c(nz) + wj1*c1(nz)) + Ft)*dt
   m%loss(CANDEP)  = sum( m%zloss(:) )

   ! Accumulate

   m%b(GRNDDEP:CANDEP) = m%b(GRNDDEP:CANDEP) + m%loss(GRNDDEP:CANDEP) 
  
   if( debug_level >1 ) then
     !if( ind == 6 ) call print_mass(ind,nSteps,nz,c)
     call print_mass(ind,nSteps,nz,c, "Inside Zdiff")
   end if

 end associate 
   
end subroutine get_mass

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module esx_ZdiffSolver

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
