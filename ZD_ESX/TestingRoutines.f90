
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!>MODULE  TestingRoutines   - analytical solutions and small help routines
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
                          module TestingRoutines
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

  !- collection of analytical solutions and small help routines for testing

!ESX  use gamma_ml, only : my_gamma   ! Rosetta-code gamma function. 
  implicit none
  private

  public :: an_sol_m1
  public :: an_sol_c0
  public :: an_sol_c1
  public :: an_sol_c2

  public :: def_prof

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! analytical solution for 1D diffusion
! m = remaining mass
! --unit ground source  !!??? No source term here
! --constant K(z)
! --deposition included
! --no upper boundary

elemental function an_sol_m1(t,K,Vd) result(m)
  real,    intent(in) :: t,K,Vd
  real :: m

  m = exp(Vd**2*t/K) * erfc(Vd*sqrt(t/K));

end function an_sol_m1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! analytical solution for 1D diffusion
! c = concentration
! --unit ground source
! --constant K(z)
! --deposition included
! --no upper boundary

elemental function an_sol_c1(t,z,K,Vd) result(c)
  real,    intent(in) :: t,z,K,Vd
  real :: c
  real :: pi

  pi = 4.0*atan(1.0)
  c = 1/sqrt(pi*K*t) * exp(-z**2/(4*K*t)) - &
      Vd/K * exp(Vd*z/K + Vd**2*t/K) * erfc(z/sqrt(4*K*t) + Vd*sqrt(t/K))
end function an_sol_c1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> Analytical solution from Hanna's Handbook on atmospheric diffusion
!! (Sect. 8.2.1), but now using a factor two to account for "reflection"
! of negative-z solutions.
! analytical solution for 1D diffusion
! c = concentration
! --instantaneous ground source of Q
! --constant K(z)
! --no upper boundary

elemental function an_sol_c0(Q,t,z,K) result(c)
  real,    intent(in) :: Q,t,z,K
  real :: c
  real :: pi
  pi = 4.0*atan(1.0)
  c =  2 * Q/sqrt(4*pi*K*t) * exp(-z**2/(4*K*t))
end function an_sol_c0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! analytical solution for 1D diffusion
! c = concentration
! -- unit ground source
! -- variable K(z) = Ka(z/za)^n
! -- no deposition
! -- no upper boundary
elemental function an_sol_c2(t,z,Ka,za,n) result(c)
  real,    intent(in) :: t,z,Ka,za, n
  real :: c
  real :: r

  r = 2-n
  c = r/(za*gamma(1.0/r)) * (za**2/(r**2*Ka*t))**(1/r) * &
    exp(-za**(2-r)*z**r/(r**2*Ka*t))
  
end function an_sol_c2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function def_prof(n,i1,i2,val) result (x)
  integer, intent(in) :: n, i1, i2
  real, intent(in) :: val
  real, dimension(n) :: x
  x(:) = 0.0
  x(i1:i2) = val
end function def_prof


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
                        end module TestingRoutines
! __________________________________________________________________________ !
! __________________________________________________________________________ !


program Testing
  use TestingRoutines
  implicit none
  integer, parameter :: NZ = 20 ! 200
  real, dimension(NZ) :: zbnd, z, dz, c0, c2
  integer :: iz, it
  real :: t, K=0.5, pi=4.0*atan(1.0), sumcdz

  !zbnd = (/ 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 32.0, 64.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0 /)
  zbnd = (/ ( 0.5 + 1.0*iz, iz=1,NZ )  /)

  z(1)  = 0.5*zbnd(1)
  dz(1) = zbnd(1)
  do iz = 1, size(zbnd)-1
    z(iz+1) = 0.5 * ( zbnd(iz+1) + zbnd(iz) ) 
    dz(iz+1) = zbnd(iz+1) - zbnd(iz) 
  end do
  
  do it = 0, 200, 20

    t=max(0.01, 1.0*it)

   ! here we run both solutions, with the same setup, using n=0 and Ka=0.5
   ! to get constant Kz=0.5 in the c2 soln,  and inst. unit emission

    c0(:) = an_sol_c0( Q=1.0, t=t, z=z(:), K=0.5 )  ! Concentrations at time t, heights z
    c2(:) = an_sol_c2(t,z(:),Ka=0.5,za=1.0,n=0.0) ! c2 with fake constant Kz=0.5

   ! Lots of extra printouts, all showing that the sum of c0*dz is 0.5, whereas
   ! the sum of c2*dz is 1.0, as it should be.

    print "(f8.2,20es12.3)", t, c0(1:8)!, dot_product(c0,dz)
  !  print "(f8.2,20es12.3)", t, c0(1:3),  c2(1:3), dot_product(c0,dz), dot_product(c2,dz)
  end do
  print "(8x,20f12.1)", z(1:8) ! , dot_product(c0,dz)
end program Testing
