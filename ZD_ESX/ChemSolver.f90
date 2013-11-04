! <ChemSolver.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2012 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                          module ChemSolver
! MOD OD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !=======================================================================!
  ! The following chemical solver uses variable chemical timesteps and
  ! is based on the scheme suggested in J.G. Verwer and D. Simpson (1995)
  ! "Explicit methods for stiff ODEs from atmospheric chemistry",
  ! Aplied Numerical Mathematics 18 (1995) 413.
  !
  ! Note that the exact formula used have been re-arranged for greater
  ! efficiency (Steffen Unger).
  ! Variable Dchem is used to keep track of changes from call to call.
  ! Note: decoupling of (NO3,N2O5), (PAN,CH3COO2), (MPAN,MACRO2)
  ! variable timestep (Peter Wind)
  !=======================================================================!

    !use ModelConstants, only: dp
    use CheckStops,      only: CheckStop
    use DefPhotolysis         ! => IDHNO3, etc.
    use ChemSpecs           ! => NSPEC_TOT, O3, NO2, etc.
!ESX    use ChemRates,   only: rct
    use Io_ml,             only : IO_LOG  !ESX, datewrite
    use ModelConstants, only: & !ESX KMAX_MID, KCHEMTOP, 
    !ESX                             dt_advec, & !ESX dt_advec_inv, &
                                 DebugCell, MasterProc, DEBUG_SOLVER,       &
                                 DEBUG_DRYRUN   !ESX   , USE_SEASALT
!ESX    use Met1d,  only:  rh
!ESX    use Chem1d, only: rcemis       !ESX was xn_old
!ESX    use Zchem, only : rcemis
    use ZchemData        !rct, rcemis for CM_Reactions
  implicit none

  private
  public  :: chemsolve        ! Runs chemical solver

     integer, parameter:: nchemMAX=1200  ! ESX 15
  integer, parameter:: NUM_INITCHEM=5    ! Number of initial time-steps with shorter dt
  !ESXreal, save::         DT_INITCHEM=20.0  ! shorter dt for initial time-steps, reduced for
  real, save::         DT_INITCHEM=1.0  ! shorter dt for initial time-steps, reduced for ESX
  integer, parameter  :: EXTRA_ITER = 1    ! Set > 1 for even more iteration
  real, parameter::  ZERO = 1.0e-30     !  Had _dp option in some versions


contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  subroutine chemsolve( k, dtstep, xn,Dchem,debug_level,itern)

    !.. In
    integer, intent(in) :: k   ! Needed for rates. ESX change?
    real, intent(in) :: dtstep  ! external time-step, e.g. 20 mins for dt_advec
    real, dimension(:), intent(inout) ::  xn    ! Nz,Nspec Concentrations, molec./cm3
    real, dimension(:), intent(inout) ::  Dchem ! Nz,Nspec Tendencies, d xn/dt due to chem
    integer, intent(in) :: debug_level          ! 0 for none, ..
    integer, intent(in), optional :: itern
    integer :: toiter

    logical, save ::  first_call = .true.

    real, parameter ::  CPINIT = ZERO ! 1.0e-30  ! small value for init

    !  Local
    integer :: ichem, iter,n        ! Loop indices
    integer, save ::  nchem, nzlev    ! No chem time-steps, and z-levels
    real    ::  dt2
    real    ::  P, L                ! Production, loss terms
    real    :: xextrapol   !help variable

    ! Concentrations : xold=old, x=current, xnew=predicted
    ! - dimensioned to have same size as "x"

    real, dimension(size(xn)) :: x, xold ,xnew   ! Working array [molecules/cm3]
    real, dimension(nchemMAX), save :: &
                        dti            &! variable timestep*(c+1)/(c+2)
                       ,coeff1,coeff2,cc ! coefficients for variable timestep

!======================================================

    if ( first_call ) then
        nzlev = size( rct, dim=2 )
        print *, "ChemSolver sizes spec x nz ", size(xn), nzlev 
       !print *, "ESX RATES ", rct(1)
       Dchem=ZERO
       call makedt(dtstep,dti,nchem,coeff1,coeff2,cc)
       if ( MasterProc ) then
           write(IO_LOG,"(a,i4)") 'Chem dts: nchemMAX: ', nchemMAX
           write(IO_LOG,"(a,i4)") 'Chem dts: nchem: ', nchem
           write(IO_LOG,"(a,i4)") 'Chem dts: NUM_INITCHEM: ', NUM_INITCHEM
           write(IO_LOG,"(a,f7.2)") 'Chem dts: DT_INITCHEM: ', DT_INITCHEM
           write(IO_LOG,"(a,i4)") 'Chem dts: EXTRA_ITER: ', EXTRA_ITER
           if(DEBUG_DRYRUN) write(*,*) "DEBUG_DRYRUN Solver"
       end if
       first_call = .false.
    endif

!======================================================


    !**  toiter gives the number of iterations used in TWOSTEP.
    !**  Use more iterations near ground:

     toiter = 3
     if ( present(itern) ) then
       toiter = toiter * itern
       print *, "ChemSolve increasing itern ", itern, toiter
     end if

   ! to get better accuracy if wanted (at CPU cost)
    
   !ESX toiter = toiter * EXTRA_ITER


    !** Establishment of initial conditions:
    !   Previous concentrations are estimated by the current
    !   minus Dchem because the current may be changed by
    !   processes outside the chemistry:

    !ESX do k = 2, KMAX_MID
!print *, "ESX TESTING ", size(x), size(Dchem(:,:,:,:),1), size(xnew), size(xn_2d,1), dti(1)*1.5
!do ichem = 1, 4
!print "(a,i3,a,10es10.2)", "ESX TESTING XN ", ichem, species(ichem)%name, xn_2d(ichem,:)
!end do
    !ESX do k =  k1, k2   !! ESX KCHEMTOP, KMAX_MID
!ESXBUG     do k =  1,  nzlev 


       xnew(:) = xn(:)

       x(:)    = xn(:) - Dchem(:)*dti(1)*1.5
!TEST
       x(:)    = max (x(:), ZERO )


       !*************************************
       !     Start of integration loop      *
       !*************************************


       do ichem = 1, nchem

          do n=1,NSPEC_TOT

             xextrapol = xnew(n) + (xnew(n)-x(n)) *cc(ichem)
             xold(n) = coeff1(ichem)*xnew(n) - coeff2(ichem)*x(n)

             xold(n) = max( xold(n), ZERO )
             x(n) = xnew(n)
             xnew(n) = xextrapol

          enddo

          dt2  =  dti(ichem) !*(1.0+cc(ichem))/(1.0+2.0*cc(ichem))

          where ( xnew(:) < CPINIT  )
             xnew(:) = CPINIT
          end where

!== Here comes all chemical reactions
!=============================================================================
          if ( DEBUG_DRYRUN ) then
            ! Skip fast chemistry
          else

            do iter = 1, toiter  !ESX (k)
!
! The chemistry is iterated several times, more close to the ground than aloft.
! For some reason, it proved faster for some compilers to include files as given below
! with the if statements, than to use loops.
!Just add some comments:
!At present the "difference" between My_FastReactions and My_SlowReactions
!is that in My_Reactions the products do not reacts chemically at all,
!and therefore do not need to be iterated.  We could have another class
!"slowreactions", which is not iterated or fewer times. This needs some
!work to draw a proper line ......


                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                   include 'CM_Reactions1.inc'
                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            end do !! End iterations

          ! slower? species

          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
           include 'CM_Reactions2.inc'
          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          end if ! DEBUG_DRYRUN

          if ( debug_level > 0 ) then !LOTS of output, careful. Useful
             do n=1,NSPEC_TOT
               if(iter==1) write(*, "(a,2i3,1x,a,10es10.2)") "ESX TESTING CHECM ", ichem,&
                  n, species(n)%name, xnew(n)
               call CheckStop( xnew(n) /= xnew(n), "ERROR!! ChemSolver:Nan found!") 
             end do !n
          end if
       end do ! ichem

       !*************************************
       !     End of integration loop        *
       !*************************************


       !**  Saves tendencies Dchem and returns the new concentrations:

            Dchem(:) = (xnew(:) - xn(:))/dtstep
if(k==5) then
  print *, "CHEM1D xxxxxxxxxxxxxxxxxxxxxxxxxx"
  do n = 1, NSPEC_TOT
    print "(a,9es15.4)", "CHEM1D "//trim(species(n)%name), xn(n), xnew(n), Dchem(n) 
  end do
end if
            xn(:) = xnew(:)

!ESXBUG    end do ! k vertical loop

  end subroutine chemsolve

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine  makedt(dtstep, dti,nchem,coeff1,coeff2,cc)

!=====================================================================
! Makes coefficients for two-step (written by Peter Wind, Febr. 2003)
! The formulas for coeff1, coeff2 and dti can be found in:
! J.G. Verwer and D. Simpson, "Explicit methods for stiff ODEs from
!  atmospheric chemistry", Aplied Numerical Mathematics 18 (1995) 413
!
! Note: It is better to take first some small steps, and then
!       larger steps, than to increase the timestep gradually.
!=====================================================================

 implicit none

 real, intent(in) :: dtstep
 real, dimension(nchemMAX),intent(out) :: dti,coeff1,coeff2,cc
 integer,                  intent(out) :: nchem

 real    :: ttot, dt(nchemMAX)
 real :: dt_init   ! time (seconds) with initially short time-steps
 integer :: i
!_________________________

  nchem=nchemMax !number of chemical timesteps inside dt_advec

   dt_init = NUM_INITCHEM*DT_INITCHEM

 ! - put special cases here:
  !NOT NEEDED? if(GRIDWIDTH_M>60000.0)nchem=15

!/ ** For smaller scales, but not tested
  if(dtstep<620.0) nchem = NUM_INITCHEM +int((dtstep- dt_init) / dt_init )

!/ Used for >21km resolution and dtstep>520 seconds:
!.. timesteps from 6 to nchem

   dt=(dtstep - dt_init )/(nchem-NUM_INITCHEM)

   dt(1:NUM_INITCHEM)=DT_INITCHEM     !.. first five timesteps

   if(dtstep<= dt_init )then
      nchem=int(dtstep/DT_INITCHEM)+1
      dt=(dtstep)/(nchem)
   endif
!/ **

   call CheckStop(dtstep<DT_INITCHEM, &
        "Error in Solver/makedt: dtstep too small!")
   call CheckStop(nchem>nchemMAX,&
        "Error in Solver/makedt: nchemMAX too small!")

   nchem=min(nchemMAX,nchem)

    if( MasterProc ) then

      !write(*,*)'Number of chemistry timesteps within one dtstep: ',nchem
      !27 format(' chem timestep ',I3,F13.6,' total: ',F13.6)

      ttot=0.0
      do i=1,nchem
         ttot=ttot+dt(i)
         !write(*,27)i,dt(i),ttot
      enddo

      !check that we are using consistent timesteps
      call CheckStop(abs(ttot-dtstep)>1.E-5, &
              "Error in Solver/makedt: dtstep and dt not compatible")

    endif

!.. Help variables from Verwer & Simpson
       cc(1)=1.0
       coeff2(1)=1.0/(cc(1)**2+2*cc(1))
       coeff1(1)=(cc(1)+1)**2*coeff2(1)
       dti(1)=((cc(1)+1)/(cc(1)+2))*dt(1)

       do i=2,nchem
         cc(i)=dt(i-1)/dt(i)
         coeff2(i)=1.0/(cc(i)**2+2.0*cc(i))
         coeff1(i)=((cc(i)+1.0)**2)*coeff2(i)
         dti(i)=((cc(i)+1.0)/(cc(i)+2.0))*dt(i)
         cc(i)=1.0/cc(i)
       enddo

end subroutine makedt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module ChemSolver
