module Functions_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves, aerodynamics (PsiM, PsiH,
! AerRes, bilinear-interpolation routines, and Polygon.
! ds and Standard Atmosphere p -> H conversion
!____________________________________________________________________
!
!** includes
!   troe - standrad chemical function
!   bilin_interpolate - generic, elemental - guessed bilinera method
!
!   Depends on: none - self-contained.
!   Language: F
!   History:
!   ds - 2000-Jan. 2001
!____________________________________________________________________
  use PhysicalConstants_ml, only : KAPPA, PI
  implicit none
  private

  public ::  troe
  public ::  Daily_cosine   ! Generates daily values of a variable
                            ! specified as a cosine curve over the year.
  public ::  Daily_sine     ! Generates daily values of a variable
                            ! specified as a sine curve over the year.
  public ::  Daily_halfsine ! Similar, but only half-sine curve (0..pi)
                            ! used. (E.g. for H2O2 in ACID versions)

  public :: Polygon         ! Used in deposition work.

  public :: StandardAtmos_kPa_2_km   ! US Standard Atmosphere conversion


 !/- Exner subroutines: ------------------------------------------------------

 public :: Exner_nd        ! (p/P0)**KAPPA
 public :: Tpot_2_T        ! Same as Exner_nd - but easier to remember
 public :: T_2_Tpot        ! Inverse as Exner_nd
 public :: Exner_tab       ! Tabulation. Must be called first


 !/- Interpolation constants

 real, private, parameter  ::    &
       PINC=1000.0              &  
      ,P0  =1.0e5               &    ! Standard pressure
      ,PBAS=-PINC

 real, save, private, dimension(131) ::  tab_exf  ! Tabulated Exner

 !/-- Micromet (Aerodynamic) routines ----------------------------------------

  public :: rh2vpd

  public :: AerRes

  public :: AerResM

  public :: PsiH

  public :: PsiM

  public :: PhiH       !ds Added 8/12/2004

  !/-- grid-handling subroutine

   public :: GridAllocate   ! allocates e.g. emissions, landuse an their
                            ! codes for each gridsquare
   character(len=30) :: errmsg  

  !/-- interpolation stuff
  public  :: bilin_interpolate                         !  "Generic" subroutine
  private :: bilin_interp_elem
  private :: bilin_interp_array

  real, public, dimension(0:1,0:1) :: wt    ! weighting factors, array version

  interface bilin_interpolate
     module procedure bilin_interp_array
     module procedure bilin_interp_elem
  end interface



  !========================================
  contains
  !========================================

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function troe(k0,kinf,Fc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! ds note - this isn't checked or optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,Fc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf	! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!	give Fc already as log(Fc)

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*Fc)

  end function troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_cosine(mean, amp, dmax, ndays) result (daily)
  !+
  !   Specifies cosine curve for a variable over a year

     real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
     integer, intent(in)  :: dmax         ! Day where maximum occurs
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: twopi                  ! Could use PhysiclConstants_ml
     twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_ml
                                          ! standalone

     do d = 1, ndays
      daily(d) = mean + amp * cos ( twopi * (d - dmax)/ ndays )
     end do

  end function Daily_cosine
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_sine(mean, amp, dmax, ndays) result (daily)
  !+
  !   Specifies sine curve for a variable over a year
  !   25/9/2002, ds, dmax redifined to be true dmax. Before it was
  !   80 and actually the day when the mean ocrrurred (spotted by hf)

     real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
     integer, intent(in)  :: dmax         ! Day where maximum occurs
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: shift                  ! Shifts sine curve to give max 
                                          ! when d = dmax
     real, save :: twopi                  ! Could use PhysiclConstants_ml
     twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_ml
                                          ! standalone
     shift = ndays/4.0

     do d = 1, ndays
      daily(d) = mean + amp * sin ( twopi * (d + shift - dmax)/ ndays )
     end do

  end function Daily_sine
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_halfsine(base, amp, ndays) result (daily)
  !+
  !   Specifies half-sine curve for a variable over a year, with
  !   values 1.0 at start and end, and max in mid-summer.

     real,    intent(in)  :: base, amp    ! Annual base and amplitude of sine
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: pi                     ! Could use PhysiclConstants_ml
     pi = 4.0 * atan(1.0)                 ! but I prefer to keep Functions_ml
                                          ! standalone

     do d = 1, ndays
      daily(d) = base + amp * sin ( pi * (ndays - d )/ ndays )
     end do

  end function Daily_halfsine
  

  !___________________________________________________________________________
  !+ subroutine which can be used with data such as emissions, landuse, where
  !  several indices area llowed per grid square
  !___________________________________________________________________________

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine GridAllocate(label,i,j,ic,ncmax,iland,&
                           ncmaxfound,land,nland, errmsg)

    !-- Checks if a country "ic" (or landuse type, lu) whose data has just 
    !   been read in has already been found within the given grid square.
    !   If not, the array "nland" is incremented by one and the
    !   country (or landuse) index added to "land".
    !
    !u7.2 
    !ds - moved LandAllocate from EmisDef since the same subrotuine can
    !     be used for landuse data also.
 
     character(len=*), intent(in) :: label   ! Type of data
     integer, intent(in) :: i,j
     integer, intent(in) :: ic        ! Index of country (lu) just read in
     integer, intent(in) :: ncmax     ! Max. no countries (lu) allowed

     integer, intent(out)   :: iland         ! Index of country in that grid
     integer, intent(inout) :: ncmaxfound    ! No. countries found so far
     integer, dimension(:,:,:), intent(inout) :: land   ! Land-codes
     integer, dimension(:,:),   intent(inout) ::nland   ! No. countries
     character(len=*), intent(out) :: errmsg   !  "ok" or not

     !ds integer :: nc, k, iland      ! local variables
     integer :: nc, k             ! local variables

       nc=nland(i,j)       ! nc = no. countries known so far
       errmsg = "ok"

       do k = 1,nc
          if( land(i,j,k) == ic) then
              iland = k        ! country is already in the list
              goto 100
          endif
       enddo

       nland(i,j) = nland(i,j) + 1    ! country is a new one
       land(i,j,nc+1) = ic
       iland=nc+1

       if( iland >  ncmaxfound) then
           ncmaxfound = iland
           write(*,*) "GridAlloc ", label, "increased ncmaxfound:",i,j,iland
           write(*,*) "GridAlloc ", label," now have:", &
                           (land(i,j,k),k=1,ncmaxfound)
           if ( ncmaxfound >  ncmax ) then
               errmsg = "GridAlloc ncmax ERROR" // label
               print *, errmsg
           endif
        endif
 100    continue

  end subroutine GridAllocate
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !___________________________________________________________________________
  !+ subroutines which can be used in 2-D interpolation
  !  - includes "generic" subroutine bilin_interpolate
  !___________________________________________________________________________

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine bilin_interp_array(xp,yp,ixp,iyp)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  real,    intent(in)  :: xp, yp     ! coordinates of point P (see fig.)
  integer, intent(out) :: ixp, iyp   ! integer coords P

  !/ Output:
  ! real, intent(out), dimension(0:1,0:1) :: wt    ! weights (see below)
  !-----------------------------------------------------------------------------
  ! This subroutine uses a bilinear interpolation method which suuplies the 
  ! weighting factors needed to estimate the value of a field at a point P 
  ! (input coords xp, yp) from the values at the nearest 4 grid points. 
  !
  ! This routine assumes that P is given in the coordinates of the field 
  ! which is being interpolated. If we define ixp = int(xp),iyp=int(yp),
  ! dx = xp - ixp, dy = yp - iyp,  we obtain a system: 
  !
  !        y'
  !        ^
  !        |
  !        0,1--------------------------1,1
  !        |                             |
  !        |                             |
  !        |                             |
  !        p1               *P(dx,dy)    p2
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        0,0 -------------------------1,0----------> x'
  !
  ! This subroutine outputs the weight to be given to the four corners
  ! using the array wt(0:1,0:1). 
  !
  ! For the bilinear interpolation we first calculate the weights associated
  ! with points p1,p2 along the y-axis, then interpolate these to P along the 
  ! x-axis
  !
  !  C(0,p1)  = (1-dy) * C(0,0)  + dy * C(0,1)
  !  C(1,p2)  = (1-dy) * C(1,0)  + dy * C(1,1)
  !  C(dx,dy) = (1-dx) * C(0,p1) + dx * C(1,p2)
  !           = (1-dx) * (1-dy) * C(0,0) +(1-dx) * dy * C(0,1) 
  !            +  dx   * (1-dy) * C(1,0) +   dx  * dy * C(1,1)
  !  i.e. Cp  
  !           = (1-dx-dy+dx.dy) * C(0,0)
  !            +(dy  -dx.dy)    * C(0,1)
  !            +(dx  -dx.dy)    * C(1,0)
  !            +(dx.dy)         * C(1,1)
  ! The "wt" array consists of the 4 coefficients of the C terms
  !
  ! Notes:
  !  - robust against P lying on either or both axis - no special cases are 
  !    needed.
  !  - assumes that field values exist at all corners. This is fine as long
  !    as we are using the method to interpolate from global fields.
  !-----------------------------------------------------------------------------
    real :: dx, dy, dxdy      ! local variables

    ixp = int(xp)
    iyp = int(yp)

    dx  = xp - ixp
    dy  = yp - iyp
    dxdy =dx * dy

    wt(0,0) = 1.0 - dx - dy + dxdy
    wt(0,1) = dy - dxdy
    wt(1,0) = dx - dxdy
    wt(1,1) = dxdy

  end subroutine bilin_interp_array

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental subroutine bilin_interp_elem(xp,yp,ixp,iyp,wt_00,wt_01,wt_10,wt_11)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  real,    intent(in)  :: xp, yp     ! coordinates of point P (see fig.)
  integer, intent(out) :: ixp, iyp   ! integer coords P
  real, intent(out)    :: wt_00, wt_01, wt_10, wt_11  ! weights, see below

  !-----------------------------------------------------------------------------
  ! method as for subroutine bilin_interp_array, but now we return scalar
  ! arguments so that the routine can be elemental. Not quite so elegant
  ! maybe, but elemental is nice.
  !  Now we have wt_00 = wt(0,0), wt_01 = wt(0,1), etc.
  ! Note the potential for error if the arguments are not called in the correct
  ! order!
  !-----------------------------------------------------------------------------
    real :: dx, dy, dxdy      ! local variables

    ixp = int(xp)
    iyp = int(yp)

    dx   = xp - ixp
    dy   = yp - iyp
    dxdy = dx * dy

    wt_00   = 1.0 - dx - dy + dxdy
    wt_01   = dy - dxdy
    wt_10   = dx - dxdy
    wt_11   = dxdy

  end subroutine bilin_interp_elem
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!d1.1    - checks for LenS, LenE = 0 introduced
!d1.3    - simplified following Margaret's suggestion: removed 
!          if-tests and possibility of EGS > 366
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!=======================================================================
function Polygon(jdayin,Ymin,Ystart,Ymax,Sday,LenS,Eday,LenE) &
result (Poly)
!=======================================================================

!     Calculates the value of a parameter Y with a polygon
!     distribution - currently LAI and g_pot

!            _____________       <- Ymax
!           /             \
!          /               \
!         /                 \
!        /                   \
!       |                     |  <- Ystart
!       |                     |
!       |                     |
!  ----------------------------- <- Ymin
!       S  S1            E1   E
!
!1.4 The following has been simplified
    

!   Inputs
    integer, intent(in) :: jdayin      !day of year
!d1.4 integer, intent(in) :: yydays    !no. days in year (365 or 366)
    real, intent(in) ::    Ymin        !minimum value of Y
    real, intent(in) ::    Ystart      !value Y at start of growing season
    real, intent(in) ::    Ymax        !maximum value of Y
    integer, intent(in) ::    Sday        !start day (e.g. of growing season)
    integer, intent(in) ::    LenS        !length of Start period (S..S1 above)
    integer, intent(in) ::    Eday        !end day (e.g. of growing season)
    integer, intent(in) ::    LenE        !length of end period (E..E1 above)

!  Output:
    real ::   Poly  ! value at day jday

! Local
    integer :: jday  ! day of year, after any co-ordinate change
    integer ::    S   !  start day
    integer ::    E   !  end day
    

    jday = jdayin
    E = Eday
    S = Sday

  ! Here we removed a lot of code associated with the leaf-age
  ! version of g_pot. 
       
    if ( jday  <  S .or. jday >  E ) then
       Poly = Ymin
       return
    end if

  !d1.3 - slightly re-written tests:

    if (jday <=  S+LenS  .and. LenS > 0 ) then

        Poly = (Ymax-Ystart) * (jday-S)/LenS  + Ystart 

    else if ( jday >=  E-LenE .and. LenE > 0.0 ) then   !d1.1 test for LenE

        Poly = (Ymax-Ystart) * (E-jday)/LenE + Ystart

    else

        Poly =Ymax
    end if
    

 end function Polygon

 !=======================================================================
 elemental function StandardAtmos_kPa_2_km(p_kPa) result (h_km)
 !=======================================================================
   implicit none

  !+ Converts pressure (kPa)  to height (km) for a US standard Atmosphere
  !  Valid up to 20 km
  !
  ! ds 27/7/2003

   real, intent(in) :: p_kPa
   real             :: h_km
   real :: t    ! Temperature (K)

   if( p_kPa > 22.632 ) then   ! = 11 km height
         ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
         !- use the power function replacament, m**n == exp(n*log m)

         t = 288.15/exp(-1.0/5.255876*log(p_kPa/101.325))
         h_km = (288.15-t)/6.5
   else
         h_km = 11.0 + log( p_kPa/22.632)/(-0.1576884)
   end if

 end function StandardAtmos_kPa_2_km

 !=======================================================================


  !--------------------------------------------------------------------
  function rh2vpd(T,rh) result (vpd_res)
  !This function is not currently in use.

    real, intent(in) ::  T    ! Temperature (K)
    real, intent(in) ::  rh   ! relative humidity (%)
    real :: vpd_res     ! vpd   = water vapour pressure deficit (Pa)

    !   Local:
    real :: vpSat       ! vpSat = saturated water vapour pressure (Pa)
    real :: arg

    arg   = 17.67 * (T-273.15)/(T-29.65)
    vpSat = 611.2 * exp(arg)
    vpd_res   = (1.0 - rh/100.0) * vpSat

  end function rh2vpd

  !--------------------------------------------------------------------
  function AerRes(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_ml.

!   Out:
    real :: Ra      ! =  aerodynamic resistance to transfer of sensible heat
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiH   = integral flux-gradient stability function for heat 
!...

    Ra = log(z2/z1) - PsiH(z2*Linv) + PsiH(z1*Linv)
    Ra = Ra/(Karman*uStar)
    Ra = min(Ra,9999.9)

  end function AerRes

  !--------------------------------------------------------------------
  function AerResM(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_ml.

!   Out:
    real :: Ra     ! =  aerodynamic resistance to transfer of momentum
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiM   = integral flux-gradient stability function for momentum
!...

    Ra = log(z2/z1) - PsiM(z2*Linv) + PsiM(z1*Linv)
    Ra = Ra/(Karman*uStar)
    Ra = min(Ra,9999.9)

  end function AerResM

  !--------------------------------------------------------------------
  function PsiH(zL) result (stab_h)
    !  PsiH = integral flux-gradient stability function for heat 
    !  Ref: Garratt, 1994, pp52-54

    ! In:
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    
    ! Out:
    real :: stab_h         !   PsiH(zL) 
    
   ! Local
   real :: x
 
    if (zL <  0) then !unstable
        x    = sqrt(1.0 - 16.0 * zL)
        stab_h = 2.0 * log( (1.0 + x)/2.0 )
    else             !stable
        stab_h = -5.0 * zL
    end if

  end function PsiH

  !--------------------------------------------------------------------
  function PsiM(zL) result (stab_m)
   !   Out:
   !   PsiM = integral flux-gradient stability function for momentum 
   !   Ref: Garratt, 1994, pp52-54

    real, intent(in) ::  zL    ! = surface layer stability parameter, (z-d)/L 
                               ! notation must be preserved         
    real :: stab_m
    real  :: x
 
    if( zL < 0) then !unstable
       x    = sqrt(sqrt(1.0 - 16.0*zL))
       stab_m = log( 0.125*(1.0+x)*(1.0+x)*(1.0+x*x) ) +  PI/2.0 - 2.0*atan(x)
    else             !stable
       stab_m = -5.0 * zL
    end if

  end function PsiM

!--------------------------------------------------------------------
  function PhiH(zL) result (phi_h)
    !  PhiH = flux-gradient stability function for heat 
    !  Ref: Garratt, 1994, pp52-54

    ! In:
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    
    ! Out:
    real :: phi_h         !   PhiH(zL) 
    
   ! Local
   real :: x
 
    if (zL <  0) then !unstable
        phi_h    = sqrt(1.0 - 16.0 * zL)
    else             !stable
        phi_h = 1 + 5.0 * zL
    end if

  end function PhiH
!--------------------------------------------------------------------

 !+
 !  Exner functions
 !
 !  Defined here as (p/P0)**KAPPA

 ! Where KAPPA = R/CP = 0.286
 ! P0 = 1.0e5 Pa

 ! CAREFUL:  The term Exner function  can also be used for CP * (p/P0)**KAPPA
 ! Hence notation Exner_nd  - non dimensional version
 !
 ! Tabulate  :
    ! defines the exner-function for every 1000 pa from zero to 1.3e+5 pa
    ! in a table for efficient interpolation (same procedure as used in
    ! the nwp-model, see mb1e.f)
 !
 ! Exner_nd returns the non-dimesnional excner function (p/p0)**R/CP
 !
 ! Added 7/4/2005, Dave, based upon tpi code from tiphys
 ! Test prog at end along with results.
 !----------------------------------------------------------------------------



  !-------------------------------------------------------------------
  subroutine Exner_tab()
  !
    real    :: p
    integer :: i

	do i = 1,131
	  p = PBAS + i*PINC
	  ! tpi(i) = CP*(p/1.0e+5)**KAPPA ! With CP!!!!
	  tab_exf(i) = (p/1.0e+5)**KAPPA  ! Without CP
	enddo

  end subroutine Exner_tab
  !-------------------------------------------------------------------

  elemental function Exner_nd(p) result(exf)

     real, intent(in) :: p    ! Pressure, p
     real :: exf, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = x1
        exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))

  end function Exner_nd
  !-------------------------------------------------------------------

  elemental function Tpot_2_T(p) result(fTpot)
     ! Identical to Exner_nd
     ! Usage:   T = Tpot * Tpot_2_T(p)

     real, intent(in) :: p    ! Pressure, p
     real :: fTpot, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = x1
        fTpot =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))

  end function Tpot_2_T
  !-------------------------------------------------------------------
  elemental function T_2_Tpot(p) result(fT)
     ! Iinvese of Exner_nd
     ! Usage:   Tpot = T * T_2_Tpot(p)

     real, intent(in) :: p    ! Pressure, p
     real :: fT, exf, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = x1
        exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))
        fT = 1.0/exf

  end function T_2_Tpot
  !-------------------------------------------------------------------

!program Test_exn
!  use Exner_ml
!  use PhysicalConstants_ml, only : KAPPA
!  implicit none
!
!  real :: p, exf1, exf2
!  integer :: i
!
!  call Exner_tab()
!
!   do i = 1, 20
!     p = 0.05 * i * 1.0e5
!     exf1 = Exner_nd(p)
!     exf2 = (p*1.0e-5)**KAPPA
!     print "(f8.3,4f12.5)", 1.0e-2*p, exf1, exf2, Tpot_2_T(p), T_2_Tpot(p)
!  end do
!end program Test_exn
    
! Results:
! p(mb)        exf1         exf2      Tpot_2_T    T_2_Tpot
!  50.000     0.42471     0.42471     0.42471     2.35455
! 100.000     0.51778     0.51778     0.51778     1.93133
! 150.000     0.58141     0.58141     0.58141     1.71997
! 200.000     0.63124     0.63124     0.63124     1.58418
! 250.000     0.67282     0.67282     0.67282     1.48629
! 300.000     0.70881     0.70881     0.70881     1.41081
! 350.000     0.74075     0.74075     0.74075     1.34999
! 400.000     0.76957     0.76957     0.76957     1.29943
! 450.000     0.79592     0.79592     0.79592     1.25641
! 500.000     0.82025     0.82025     0.82025     1.21913
! 550.000     0.84291     0.84291     0.84291     1.18637
! 600.000     0.86414     0.86414     0.86414     1.15722
! 650.000     0.88414     0.88414     0.88414     1.13105
! 700.000     0.90307     0.90307     0.90307     1.10734
! 750.000     0.92105     0.92105     0.92105     1.08571
! 800.000     0.93820     0.93820     0.93820     1.06587
! 850.000     0.95461     0.95461     0.95461     1.04755
! 900.000     0.97033     0.97033     0.97033     1.03058
! 950.000     0.98544     0.98544     0.98544     1.01477
!1000.000     1.00000     1.00000     1.00000     1.00000
end module Functions_ml
