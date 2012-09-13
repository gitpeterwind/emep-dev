
! <Solver.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
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

                          module Chemsolver_ml
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

    use Aqueous_ml,        only: aqrck, ICLOHSO2, ICLRC1, ICLRC2, ICLRC3
    use CheckStop_ml,      only: CheckStop
    use DefPhotolysis_ml         ! => IDHNO3, etc.
    !ESX use EmisDef_ml,        only: QSSFI, QSSCO, QDUFI, QDUCO, QPOL, &
    !ESX                              QROADDUST_FI, QROADDUST_CO
    use Emissions_ml,      only: KEMISTOP
    use ChemGroups_ml,     only: RO2_POOL, RO2_GROUP
    use ChemSpecs_tot_ml           ! => NSPEC_TOT, O3, NO2, etc.
    use Chemfields_ml, only : NSPEC_BGN  ! => IXBGN_  indices and xn_2d_bgn
    use ChemRates_rct_ml,   only: rct
    !ESX use ChemRates_rcmisc_ml,only: rcmisc
    use GridValues_ml,     only : GRIDWIDTH_M
    use Io_ml,             only : IO_LOG, datewrite
    use ModelConstants_ml, only: KMAX_MID, KCHEMTOP, dt_advec,dt_advec_inv, &
                                 DebugCell, MasterProc, DEBUG_SOLVER,       &
                                 DEBUG_DRYRUN, USE_SEASALT
    use Par_ml,            only: me, MAXLIMAX, MAXLJMAX
    use PhysicalConstants_ml, only:  RGAS_J
    use Setup_1dfields_ml, only: rcemis,        & ! photolysis, emissions
                                 xn_2d,         &
                                 rh,            &
                                 Fgas,   & ! fraction in gas-phase, for SOA
                                 amk
                                 !FUTURE rcnh3,         & ! NH3emis
 use Setup_1dfields_ml,     only : itemp, tinv, rh, x=> xn_2d, amk
    use ChemFunctions_ml, only :VOLFACSO4,VOLFACNO3,VOLFACNH4 !TEST TTTT
  implicit none

  private
  public  :: chemistry        ! Runs chemical solver

  INCLUDE 'mpif.h'

  integer::  STATUS(MPI_STATUS_SIZE),INFO
     integer, parameter:: nchemMAX=15
  integer, parameter:: NUM_INITCHEM=5    ! Number of initial time-steps with shorter dt
  real, save::         DT_INITCHEM=20.0  ! shorter dt for initial time-steps, reduced for
  integer, parameter  :: EXTRA_ITER = 1    ! Set > 1 for even more iteration


contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  subroutine chemistry(i,j,debug_flag)

    !.. In
    integer, intent(in) ::  i,j       ! Coordinates (needed for Dchem)
    logical, intent(in) :: debug_flag

    real, dimension(:,:,:,:), save,allocatable :: &
                    Dchem  ! Concentration increments due to chemistry

    logical, save ::  first_call = .true.

    real, parameter ::  CPINIT = 0.0 ! 1.0e-30  ! small value for init

    !  Local
    integer, dimension(KCHEMTOP:KMAX_MID) :: toiter
    integer ::  k, ichem, iter,n    ! Loop indices
    integer, save ::  nchem         ! No chem time-steps
    real    ::  dt2
    real    ::  P, L                ! Production, loss terms
    real    :: xextrapol   !help variable

    ! Concentrations : xold=old, x=current, xnew=predicted
    ! - dimensioned to have same size as "x"

    real, dimension(NSPEC_TOT)      :: &
                        x, xold ,xnew   ! Working array [molecules/cm3]
    real, dimension(nchemMAX), save :: &
                        dti             ! variable timestep*(c+1)/(c+2)
    real, dimension(nchemMAX), save :: &
                        coeff1,coeff2,cc ! coefficients for variable timestep

!======================================================

    if ( first_call ) then
       allocate( Dchem(NSPEC_TOT,KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX))
       Dchem=0.0
       call makedt(dti,nchem,coeff1,coeff2,cc)
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

    toiter(KCHEMTOP:5)        = 1    ! Upper levels - slow chemistry
    toiter(6:KEMISTOP-1)      = 2    ! Medium and cloud levels
    toiter(KEMISTOP:KMAX_MID) = 3    ! Near-ground, emis levels

   ! to get better accuracy if wanted (at CPU cost)
    toiter = toiter * EXTRA_ITER


    !** Establishment of initial conditions:
    !   Previous concentrations are estimated by the current
    !   minus Dchem because the current may be changed by
    !   processes outside the chemistry:

    do k = 2, KMAX_MID

       xnew(:) = xn_2d(:,k)

       x(:)    = xn_2d(:,k) - Dchem(:,k,i,j)*dti(1)*1.5
       x(:)    = max (x(:), 0.0)


       !*************************************
       !     Start of integration loop      *
       !*************************************


       do ichem = 1, nchem

          do n=1,NSPEC_TOT

             xextrapol = xnew(n) + (xnew(n)-x(n)) *cc(ichem)
             xold(n) = coeff1(ichem)*xnew(n) - coeff2(ichem)*x(n)
             xold(n) = max( xold(n), 0.0 )
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

            do iter = 1, toiter(k)
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

                !if(k>=KCHEMTOP)then

                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                   include 'CM_Reactions1.inc'
                   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

                !endif
                !if(k>=6)then
                !   include 'My_FastReactions.inc'
                !endif
                !if(k>=KEMISTOP)then
                !   include 'My_FastReactions.inc'
                !endif
            end do !! End iterations
          ! Just before SO4, look after slower? species

          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
           include 'CM_Reactions2.inc'
          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          end if ! DEBUG_DRYRUN

       end do ! ichem

       !*************************************
       !     End of integration loop        *
       !*************************************


       !**  Saves tendencies Dchem and returns the new concentrations:

            Dchem(:,k,i,j) = (xnew(:) - xn_2d(:,k))*dt_advec_inv
            xn_2d(:,k) = xnew(:)


    enddo ! End of vertical k-loop

  end subroutine chemistry

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine  makedt(dti,nchem,coeff1,coeff2,cc)

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
  if(dt_advec<620.0) nchem = NUM_INITCHEM +int((dt_advec- dt_init) / dt_init )

!/ Used for >21km resolution and dt_advec>520 seconds:
!.. timesteps from 6 to nchem

   dt=(dt_advec - dt_init )/(nchem-NUM_INITCHEM)

   dt(1:NUM_INITCHEM)=DT_INITCHEM     !.. first five timesteps

   if(dt_advec<= dt_init )then
      nchem=int(dt_advec/DT_INITCHEM)+1
      dt=(dt_advec)/(nchem)
   endif
!/ **

   call CheckStop(dt_advec<DT_INITCHEM, &
        "Error in Solver/makedt: dt_advec too small!")
   call CheckStop(nchem>nchemMAX,&
        "Error in Solver/makedt: nchemMAX too small!")

   nchem=min(nchemMAX,nchem)

    if( MasterProc ) then

      write(*,*)'Number of chemistry timesteps within one dt_advec: ',nchem
      27 format(' chem timestep ',I3,F13.6,' total: ',F13.6)

      ttot=0.0
      do i=1,nchem
         ttot=ttot+dt(i)
         write(*,27)i,dt(i),ttot
      enddo

      !check that we are using consistent timesteps
      call CheckStop(abs(ttot-dt_advec)>1.E-5, &
              "Error in Solver/makedt: dt_advec and dt not compatible")

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

end module Chemsolver_ml
