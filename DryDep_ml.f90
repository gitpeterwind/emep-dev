
module DryDep_ml

  ! Module contains relevant parts of old readpar_mach and drydep.f
  ! MACHO method. This is the drag-coefficient based approach of
  ! BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

 !-- model specific dry dep values are set in My_DryDep_ml:

 use My_DryDep_ml, only : Init_vd, & ! Initialisation
                          SEA, &     ! Value of ind for sea squares
                          NDRYDEP_CALC, &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          vd_min, &  ! Min. Vd value
                          vd_day, &  ! Additional daytime Vd
                          DepLoss, Add_ddep, &
                          Dep        ! Mapping (type = depmap)

 use Chemfields_ml ,       only : xn_adv,cfac
 use GenSpec_adv_ml,       only : NSPEC_ADV
 use GridValues_ml ,       only : GRIDWIDTH_M,xmd,xm2,carea,iclass, &
                                  rpol2,       &  !u3 for seasonal  Vg
                                  i_glob, j_glob   !u1 for testing
 use MassBudget_ml,        only : totddep
 use Met_ml,               only : roa,fm,fh,z_bnd,th2m,th,ps,u,v,z_mid&
                                  ,ustar,foundustar
 use ModelConstants_ml,    only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  daynumber,  &  !u3 for seasonal  Vg
                                  ATWAIR, atwS, atwN
 use Par_ml,               only : me, li0,li1,lj0,lj1
 use PhysicalConstants_ml, only : XKAP, PI
 use Radiation_ml,         only : zen       ! zenith angle (degrees)
 implicit none
 private

 public :: drydep

  logical, private, save :: my_first_call = .true.
  logical, private, parameter :: DEBUG_VG = .false.
 !integer, private, parameter :: debug_i=79, debug_j=56 ! Eskdalemuir
 integer, private, parameter :: debug_i=97, debug_j=62 !  Waldhof
 !integer, private, parameter :: debug_i=73, debug_j=48 ! Mace Head
 !integer, private, parameter :: debug_i=91, debug_j=71 ! Rorvik


 contains
  subroutine drydep

 real, dimension(NDRYDEP_CALC) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 1m
         ,vg_fac         ! Loss factor due to dry dep.

 integer i, j, n, ncalc, nadv, ind   ! help indexes

 real veff,     & ! Effectice dry dep. at the centre of layer
      vd1m,     & ! Dry dep. at 1 m
      xnxh,     & ! argument in the exp. decay of dry dep.
      xnxi        ! amount dry deposited

 real ustar_sq, & ! u star squared
      abshd,    & ! abs value of surfae flux of sens. heat W/m^2
      uvhs,     & ! horizontal vind
      inv_CdV     ! 1/(Cd * V) where Cd is drag coefficient and V is wind-speed
                  ! for downscaling dry dep. according to atm. conditions
                  ! was reduction_factor

 real convfac,  & ! rescaling to different units
      convfac2, & ! rescaling to different units
      dtz,      & ! scaling factor for veff ( = dt/z, where dt=timestep and 
                  ! z = height of layer)
      dayfac      ! Used to distinguish day/night deposition rates

  real :: sin2, cos2, seasonfac   !u3 For seasonal variation

  real, save :: inv_gridarea  ! inverse of grid area, m2




!     first calculate the 1m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

  if ( my_first_call ) then 

     call Init_vd()
     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

     my_first_call = .false.

  end if !  my_first_call

  ustar_sq = -99.99   ! for debugging

 !u3 ========= time facs
 !u3  - seasonal variation used, taken from Lagrangian O3 model which
 !      was taken from the Lag. Nox model.....

 !.. Variation of deposition vel. with fi and time of year.
        sin2 = sin((daynumber-32.0)*PI/365.25)
        sin2 = sin2*sin2
        cos2 = 1 - sin2

 !========= time facs

  do j = lj0,lj1
    do i = li0,li1

      ind = iclass(i,j)     ! land-use index

!eb
!eb     Calculating the drag-coefficient and inv_CdV (=1/Cd.V)
!eb
      abshd = abs(fh(i,j,1))       ! Heat flux

     !-- we calculate a reudtcion_factor, which is equal to the drag
     !   coefficient (C_H in BJ98) x wind speed (v_H)

      if (abshd < 0.0999999) then
         if(foundustar)then                      !pw u3
            ustar_sq = ustar(i,j,1)*ustar(i,j,1)
         else
          ustar_sq = fm(i,j,1)/roa(i,j,KMAX_MID,1)
         endif
          ustar_sq = max(ustar_sq, 1.e-10)
          uvhs = 0.5*sqrt( (u(i,j,KMAX_MID,1)+u(i-1,j,KMAX_MID,1))**2  &
                     + (v(i,j,KMAX_MID,1)+v(i,j-1,KMAX_MID,1))**2 )
          inv_CdV = 0.74*uvhs/ustar_sq
      else
          inv_CdV = (th(i,j,KMAX_MID,1)/th2m(i,j,1) - 1.) &
                      *ps(i,j,1)/(XKAP*abshd)
          inv_CdV = abs(inv_CdV)
      endif

    !u1
    ! We introduce a cap on inv_CdV after finding that sometimes it
    ! can have values of several thousands. This leads to a very
    ! great gradient even for gases with almost no deposition (e.g.
    ! ozone over oceans). I (ds) don't trust HIRLAM enough to allow that....

      inv_CdV = min(inv_CdV,500.0)

      convfac = (ps(i,j,1) - PT)*carea(KMAX_MID)*xmd(i,j)/ATWAIR

      !BUG??  dtz      = dt_advec/z_mid(i,j,KMAX_MID)

      dtz      = dt_advec/z_bnd(i,j,KMAX_BND-1)

      if ( zen(i,j) < 90 ) then
          !u3 dayfac = 0.8    !! Approx. of daytime average
          !u6 dayfac = 1.0    !! as in Hov et al, to get max Vd
          dayfac = 0.8    !!  Approx. of daytime average
      else
          !y3 dayfac = 0.0
          !dayfac = 0.25    ! stomatal closure as in Lagr. model
          dayfac = 0.0    ! stomatal closure as in Lagr. model
      end if

      ! in GridValues rpol2 is the square of (r/AN), so we need the sqrt:

      !u6b seasonfac = dayfac * sin2 + sqrt(rpol2(i,j)) * cos2
      seasonfac = sin2 + sqrt(rpol2(i,j)) * cos2




      if ( DEBUG_VG .and. i_glob(i)==debug_i .and. j_glob(j)==debug_j) then
          print "(a30,4i4)", "DEBUG DryDep me, i,j, ind ", me, i,j, ind
          print "(a20,2f8.2)", "DEBUG zen, dayfac", zen(i,j), dayfac
          print "(a20,2f8.2)", "DEBUG sin2, cos2",  sin2, cos2
          print "(a20,2f8.2)", "DEBUG rpol,seasonfac",  seasonfac,sqrt(rpol2(i,j))
         !print "(a24,4f8.2)", "DEBUG dt, z, dtz, red", &
         !       dt_advec, z_bnd(i,j,KMAX_BND-1), dtz, inv_CdV
         !print "(a20,es15.5)", "DEBUG u*^2: ", ustar_sq
         !print "(a20,3f8.2,es12.4,f8.3)", "DEBUG met: ", &
         !   fh(i,j,1), th(i,j,KMAX_MID,1),th2m(i,j,1), &
         !        fm(i,j,1),roa(i,j,KMAX_MID,1)
      end if


      do ncalc = 1, NDRYDEP_CALC
!
!     integration of the dry dep.removal from the lowest model layer.
!                    
        !u3 ========= time facs
        !u3 - until a new dep module is ready, we return to the early EMEP
        !u3 method  proposed by  Lemhaus et al., Hov et al., 1988
        !   (used in the EMEP Lagr. O3 model, surprise surprise...)

! hf correction July 1

         if ( IND == SEA ) then
            vd1m = vd_min(ncalc,ind) + vd_day(ncalc,ind)
         else
            vd1m = vd_min(ncalc,ind) + dayfac * seasonfac * vd_day(ncalc,ind)
         endif


        !u3 vd1m = vd_min(ncalc,ind) + dayfac * vd_day(ncalc,ind)
        !u6b test:

        !u3 just make 1 dep vel for now, as in Hov et al. :

        !u6b: vd1m = vd_min(ncalc,ind) + vd_day(ncalc,ind)

       !u3 ========= time facs over land squares

       !u3 ========= time facs

        gradient_fac(ncalc) = 1./(1. + vd1m*inv_CdV)
        veff = vd1m*gradient_fac(ncalc)

        vg_fac (ncalc) = 1.0 - exp ( -veff * dtz ) 

        ! Rorvik specials....
        if ( DEBUG_VG .and. & ! ncalc == 5 .and. &   ! 5 is NO2
                 i_glob(i)==debug_i .and. j_glob(j)==debug_j) then
          print "(a15,i4,6f10.4)", "DEBUG veff ", ncalc, &
              100.0*vd1m, sqrt(rpol2(i,j)), seasonfac, &
                   gradient_fac(ncalc),  100.0*veff, vg_fac(ncalc)
        end if
            
    end do ! n


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

      do n = 1, NDRYDEP_ADV 
         nadv    = Dep(n)%adv
         ncalc   = Dep(n)%calc

         DepLoss(nadv) =   vg_fac( ncalc )  * xn_adv( nadv,i,j,KMAX_MID)

         xn_adv( nadv,i,j,KMAX_MID) = &
             xn_adv( nadv,i,j,KMAX_MID) - DepLoss(nadv)

         cfac(nadv, i,j) = gradient_fac( ncalc )

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac

         if ( DEBUG_VG .and. i_glob(i)==debug_i .and. j_glob(j)==debug_j) then
              print "(a30,3i4)", "DEBUG DryDep n, adv, calc ",  n,nadv, ncalc
              print "(a20,f12.5)", "DEBUG grad fac", gradient_fac( ncalc)
              print "(a20,2e12.4)", &
                "DEBUG xn, DepLoss ", xn_adv(nadv, i,j,KMAX_MID), DepLoss(nadv)
              print "(a20,2f8.4)", "DEBUG gv_fac( ncalc)", &
                 vg_fac(ncalc), 1.0-vg_fac(ncalc)
        end if

       end do ! n

       ! inv_gridarea = xm2(i,j)/(GRIDWIDTH_M*GRIDWIDTH_M)
       convfac2 = convfac * xm2(i,j) * inv_gridarea


      !.. Add DepLoss to budgets if needed:

       call Add_ddep(i,j,convfac2)

     enddo
   enddo
 end subroutine drydep

end module DryDep_ml
