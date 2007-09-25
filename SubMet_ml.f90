module SubMet_ml
!=============================================================================
!+
! Description
!  Module for reading NWP data from file, and for calculating sub-grid
!  meteorology for each land-use. This module is designed to mimic to
!  some extent modules from the full EMEP model, but is not fully compatible.
!  The sub-grid part of this module is also undergoing constant change!!
!=============================================================================

!ToDo: do something with ios

use MicroMet_ml, only :  PsiH, PsiM, AerRes    !functions
use Io_ml, only : IO_STAB, open_file
use PhysicalConstants_ml, only : PI, RGAS_KG, CP, GRAV, KARMAN, CHARNOCK, T0
use TimeDate_ml, only: current_date

implicit none
private

!Subroutines

public :: Get_Submet    ! calculates met. data for sub-grid areas


logical, private, parameter ::  DEBUG_SUB = .false.  ! for extra tests/printouts

contains
!=======================================================================
!
!  subroutine Get_NWPmet()
!
!!---------------------------------------------------------------
!!   Reads in one time step from metdata file, calculates  stability
!!   factors, etc.
!!---------------------------------------------------------------
!
!
!  !..Local
!    integer :: ios  ! iostat variable
!   
!  ! Sigma co-ordinate of lowest layer top
!    real, parameter :: SIGMA = 0.994   ! from vgrid2.out and courant96 prog.
!    
!    real, parameter :: P100=1.0e4   ! pressure at domain-top (100 hPa) in Pa
!    real, parameter :: PT = 10000.0 ! PT = 10000.0  (100 mb) in Pa
!
!    real :: pp, T_ref, rho_ref
!
!
!!   ******************
!!    call read_metdata(ios)
!
!    psurf = psurf * 100.0   ! hPa -> N/m^2
!    Hd = -Hd           ! normal convention now, Hd away from surface!
!    Hl = -Hl
!    cl     = cl     * 0.01            ! convert percent to fractions
!    cllow  = cllow  * 0.01
!    clmed  = clmed  * 0.01
!    clhigh = clhigh * 0.01
!    theta_2m = t2 / (psurf*1.0e-5)**XKAP    ! get pot. temp. (in K) at 2m
!                                          
!!   ******************
!                                           
!
!!   Here we calculate parameters at the centre of the grid-cell 
!!   (ca. 45m), which we use as our reference height.
!
!
!
!      pp      = PT + SIGMA*(psurf-PT)                ! pressure at mid-cell
!    
!      T_ref   = theta_ref*exp(XKAP*log(pp*1.0e-5))    ! temp. at mid-cell
!
!      rho_ref = pp/(RGAS_KG*T_ref)                       ! density at mid-cell
!      z_ref   = (psurf-pp)/(rho_ref*GRAV)       ! height of mid-cell
!
!
!    !   we must use L (the Monin-Obukhov length) to calculate deposition,
!    !   therefore we calculate u*, t* from NWP-model data. 
!
!
!      rho_surf       = psurf/(RGAS_KG*t2)      ! at surface
!
!    !ds - new assumption - get u_ref directly from NWP model!
!
!      u_ref = u
!
!      ustar_nwp = sqrt(tau/rho_surf)
!      ustar_nwp = max(ustar_nwp,1.0e-2)  ! prevents zero values of ustar_nwp
!
!      tstar_nwp = -Hd/ ( CP*rho_surf*ustar_nwp )
!      
!      invL_nwp  = KARMAN * GRAV * tstar_nwp/(ustar_nwp * ustar_nwp * theta_2m)
!
!      !.. we limit the range of 1/L to prevent numerical and printout problems
!      !   This range is very wide anyway.
!
!      invL_nwp  = max( -1.0, invL_nwp ) !! limit very unstable
!      invL_nwp  = min(  1.0, invL_nwp ) !! limit very stable
!
!    !Note that in the above equation the NWP estimate of temperature at 2m
!    !above the ground is taken as an approximation for surface temperature.
!
!  end subroutine Get_NWPmet
!
!=======================================================================
  subroutine Get_Submet(h,t2,Hd,LE,psurf, z_ref, u_ref, qw_ref, & ! in-
                        debug_flag,                        &    ! in
                        ustar, invL,                       &    ! in-out
                        z0,d, Ra_ref,Ra_3m,rh,vpd)                    ! out

!---------------------------------------------------------------
!  Sub-grid calculations of  stability, Ra and ustar for this landuse
!---------------------------------------------------------------
!
!..The profile manipulation is introduced to calculate different Ra and Rb
!..values for different landuse types (e.g. from SEI). Thus, different z0,
!..are assumed within each EMEP square, from which one averaged value is 
!..available as the basic input from the NWP-model. Therefore, we first have a
!..grid square average for u*, T*, L and z0 (i.e. from the NWP-model).  Then,
!..we calculate the wind speed u from the wind profile (according to the
!..M-O similarity) and the NWP-model data, and from this u we calculate
!..new u* values for each z0, i.e. landuse, within each grid square. Here,
!..we have to use the old L, even although L=f(u*), since the new u* is
!..unknown before it is calculated. (An iteration process is possible.) 
!..Anyway, by using the new u*'s, we obtain
!..new T*'s, L's for each landuse type and finally Ra's and ustars for each
!..subgrid area.
!
!.. For further details, see EMEP report 3/95 ...
!
! **sc: The subroutine Get_Submet also generates output for 
! two terms utilized in Rsurface, namely rh (the relative humidity
! term required for the evaluation of the stomatal compensation point 
! point for NH_3) and vpd (the vapour pressure deficit term required
! for the evaluation of stomatal conductance).

!   Abbreviation in documentation: "tmp" abbreviates "temporary"
!----------------------------------------------------------------- 
!..In

   real, intent(in) :: h           ! height of vegetation (m)
                                   ! (set negative for sea areas)
   real, intent(in) :: t2          ! Surface (2m) temp in deg.K
   real, intent(in) :: Hd          ! Sensible Heat flux   ! Check SIGN!!!
   real, intent(in) :: LE          ! Latent   Heat flux   ! Check SIGN!!!
   real, intent(in) :: psurf       ! Pascal
   real, intent(in) :: z_ref       ! Height of mid-cell used as ref ht, ca. 45m
   real, intent(in) :: u_ref       ! Wind speed at ref ht 
   real, intent(in) :: qw_ref      ! Spec. humidity at ref. ht
   logical, intent(in) :: debug_flag   ! set true for wanted grid square

! In-Out
    real, intent(inout) :: ustar            ! u* for sub-grid   (m/s)
    real, intent(inout) :: invL             ! inverse of Monin-Obukhov length

! Out
 
    real, intent(out) :: z0               ! roughness length (m)
    real, intent(out) :: d                ! zero-plane displacament height (m)
    real, intent(out) :: Ra_ref           ! resistances for SEI  landuse
    real, intent(out) :: Ra_3m            ! resistances for SEI  landuse
    real, intent(out) :: rh               ! relative humidity
    real, intent(out) :: vpd              ! vapour pressure deficit in leaf    

!.. Local
    real :: z_refd                 ! z_ref - d   (m)
    real :: rho_surf               ! Density at surface (2 m), kg/m3
    real :: tstar                  ! T* for sub-grid,  pot. temp for 2 m
    real :: ustar_Nwp, invL_nwp    !  Store NWP values for testing
!u7.lu    real :: d                      ! displacement height ~ 0.7 * h 
    real :: Psi_h0, Psi_h2         ! stability functions, heat
    real :: tmpinvL, tmpinvL_nwp   ! 1/L values for printout
    real :: z_1m                   !  1m above vegetation
    real :: z_3m                   !  3m above ground, or top of trees
    real :: z_3md                  !  minus displacemt ht.
    
    
    logical, save ::  my_first_call = .true.
    integer, parameter ::  NITER = 1           ! no. iterations to be performed

    integer :: iter,ios                 ! iteration variable and iostat variable                     



   ! For vapour pressure calculations

    real, parameter :: ESAT0=611.0   ! saturation vapour pressure at 
                                     ! T=0 deg. C (Pa)
      
    real :: qw                       ! specific humidity (kg/kg) corrected  
                                     ! down to z_0+d metres above the ground 
    real :: esat    ! saturation vapour pressure  (Pa)
    real :: e       ! vapour pressure at surface
    real :: olde, oldvpd !.. tmp for printout
    real :: Ra_2m   ! to get 2m qw

    ! initial guesses for u*, t*, 1/L

    ustar_nwp = ustar   ! Save first values
    invL_nwp  = invL    ! Save first values

!hj .. the zero-plane displacement (d) is the height that
!     needs to be added in order to make the profile theories
!     work over a tall vegetation (see e.g. Stull (1988), p.381). 
!     This corresponds to moving our co-ordinate system upwards 
!     by a distance d.
! 
!     by definition: u(z0+d) = 0
!
!     it has been observed that d is approximately 0.7 times
!     the mean height of the vegetation (h) and z0=h/10
!     (see e.g. Sutton (1990), Ph.D. thesis, p.46).
!     The reference height for u* transformation is then
!     taken arbitrarily at about 45m, the height of the
!     centre of the EMEP grid cell.
          

!.. for the landuse type water (has h=-99), we introduce a new zero plane 
!   displacement and use the Charnock relation to calculate the z0-values
! Jun 2002. ds - addded max 1cm limit for z0 over sea, because of problems
! caused by z0>1m. Garratt (section 4.1, Fig 4.2) suggested that Charnock's
! relation is only valid for  u* < 1 m/s, which gives z0 < 1 cm/s.


        if ( h  < 0.0 ) then ! water
             d  = 0.0
             z0 = CHARNOCK * ustar * ustar/GRAV
             z0 = max(z0,0.01)  ! u7.5vgc 
             z_1m   = 1.0       ! 1m above sea surface
             z_3m   = 3.0       ! 3m above sea surface
        else
             d  =  0.7 * h
             z0 =  0.1 * h
             z0 = max(z0,0.001) !  Fix for deserts, ice, snow (where, if the 
                                !  ground is bare, h=0 and hence z0=0)
             z_1m   = (h + 1.0) - d !  1m above vegetation, rel to displace ht.
             z_3m   = max(3.0,h)    !  3m above ground, or top of trees 
        end if
          
        z_refd = z_ref - d          !  minus displacement height
        z_3md  = z_3m  - d          !  minus displacement height


    do iter = 1, NITER 

        ! ****
        !   PsiM calculates the stability functions for momentum
        !   at heights z_ref (about 45m) & z0
        ! ****               
        !..calculate friction velocity based first on NWP-model PsiM-values 
        !..and u_ref. The NWP-model PsiM-values are used despite the fact that
        !..L=F(u*), since we do not know the EMEP subgrid averaged 
        !..z0-values ...

        if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
            write(6,"(a12,i2,5f8.3,f12.3)") "UKDEP SUBI", iter, &
                       h, z0, d, z_refd, z_3md, invL
        end if

       ustar = u_ref * KARMAN/ &
         (log( z_refd/z0 )-  PsiM( z_refd*invL ) + PsiM( z0*invL ) )

       ustar = max(ustar,1.0e-2)

    !  We must use L (the Monin-Obukhov length) to calculate deposition,
    ! Thus, we calculate T* and then L, based on sub-grid data. 


        rho_surf = psurf/(RGAS_KG * t2 )


    ! New 1/L value ....

        !u7.lu invL =  KARMAN * GRAV * tstar/(ustar * ustar * theta_2m)
        !u7.lu - don't worry about diff between pot temp and tempf or now:

        invL =  -KARMAN * GRAV * Hd / ( CP*rho_surf*ustar*ustar*ustar * t2)

      !.. we limit the range of 1/L to prevent numerical and printout problems
      !   This range is very wide anyway.

        invL  = max( -1.0, invL ) !! limit very unstable
        invL  = min(  1.0, invL ) !! limit very stable

    end do ! iter


    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        write(6,"(a12,4f8.3)") "UKDEP SUBL", z0, d, z_refd, invL
    end if

    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
       !!  modulo( met_dd, PRINT_DD)  == 0 ) then  ! print out every 10 days
         if ( my_first_call ) then ! title line

                call open_file(IO_STAB,"w","Stab.out",needed=.true.)

                write(unit=IO_STAB, fmt="(3a3, a6, 3a8,2a7, 2a6)") &
                 "mm", "dd", "hh", "t2_C", "Hd", &
                 "L_nwp", "L  ", "z/L_nwp", "z/L ", "u*_nwp", "u*"
                my_first_call = .false.
         end if

            !/ Limit 1/L values for printout
            !tmpinvL_nwp = min( 1.0/invL_nwp, 9999.9)
            !tmpinvL_nwp = max( tmpinvL_nwp, -9999.9)
            !tmpinvL     = min( 1.0/invL    , 9999.9)
            !tmpinvL     = max( tmpinvL    , -9999.9)

            write(unit=IO_STAB, &
              fmt="(3i3, f6.1, 3f8.2, 2f7.2, 2f6.2)") &
              current_date%month, &
              current_date%day, &
              current_date%hour, &
              t2-273.15, Hd, invL_nwp, invL, &
              z_refd*invL_nwp, z_refd*invL, ustar_nwp, ustar
    end if



!     *** Aerodynamic resistances for each landuse type k ***
!      Ra_ref is used to estimate the aerodynamic resistance to latent
!      heat transfer from height z_ref to z0+d and from height  h+3 to 
!      z0+d, respectively.
!      Only Ra_ref is used in the present subroutine.  
      
        Ra_ref = AerRes(z0,z_refd,ustar,invL,KARMAN)
        Ra_3m  = AerRes(z0,z_3md,ustar,invL,KARMAN)
        Ra_2m  = AerRes(z0,1.0+z_1m,ustar,invL,KARMAN)

    if ( Ra_ref < 0 .or. Ra_3m < 0 .or. Ra_2m < 0  ) then
      print *, "RAREF NEG ", z0, z_refd, ustar, invL, KARMAN
    end if
    if ( DEBUG_SUB .and.  Ra_3m > Ra_ref ) then
        print "(a22,f12.3)", "ERROR!!! Ra_ref<Ra_3",  Ra_ref, Ra_3m
    end if
    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        write(6,"(a22,f12.3)") "UKDEP SUB5 Ra_ref",  Ra_ref
    end if


!  *****  Calculate rh and vpd  *********

!  
!....The model has the specific humidity qw_ref for the lowest model layer as 
!    an input obtained from the NWP model.  We now correct this down to 
!    z_0+d metres above the ground using q(z_0+d) = q(z2) + Q*Ra_ref,
!    where Ra_ref is the same as calculated above. Q is latent heat flux 
!    (W/m^2) divided by latent heat of vap. (2.5e6  J/kg). From input data, 
!    Hl is directed downwards (!), so  we used Q = -Hl/2.5e6
!    but above ds reversed it, so we go back to the normal micromet. 
!    formula Q = Hl/2.5e6

                       
!FIX        qw = qw_ref + Hl/2.5e6 * Ra_ref          ! careful with signs!
        qw = qw_ref  + LE/2.5e6 * ( Ra_ref - Ra_2m)  ! 2m qw


      !..   qw is in kg/kg  so  e = qw*psurf/epsilon
      !..   to get e in Pascal.

       e = qw * psurf/0.622

!The equation below relies on the following assumptions: epsilon=0.622
!and L=2.5x10^6 Joules/Kg, where "epsilon" and "L" denote the ratio of 
!the gas constants for dry air against water vapour and the latent
!heat of vaporization, respectively.


     esat = ESAT0 * exp(0.622*2.5e6*((1.0/T0) - (1.0/t2))/RGAS_KG )

    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        print "(a15,2f12.6,2f12.3)", "UKDEP SUB water", qw_ref, qw, LE, 100.0*e/esat
    end if

   !ds - fix for bug found by Svetlana:.
   ! Straighforward calculation sometimes gives rh<0 or rh>1.0 -
   ! probbaly due to mismatches between the assumptions used for the stability
   ! profile here and in HIRLAM. Here we set crude limits on e to prevent
   ! impossible rh values at least:


     e = max(0.001*esat,e)    ! keeps rh >= 0.1%
     e = min(esat,e)          ! keeps rh <= 1
     rh = e/esat

!ds  ****  leaf sat. vapour pressure

      vpd    =  0.001*(esat-e)     ! gives vpd in kPa !
      vpd    =  max(vpd,0.0) 


      !if ( MY_DEBUG ) then             ! useful for comparisons
      !  olde   = qw_ref * psurf/0.622
      !  oldvpd = 0.001*(esat-olde)
      !  oldvpd =  max(oldvpd,0.0) 
      !end if
      !esat =  0.tab_esat_Pa( t_ref )
    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        write(6,"(a22,2f12.4)") "UKDEP SUB7 e/esat, rh", e/esat, rh
    end if

  end subroutine Get_Submet
! =====================================================================

end module SubMet_ml
