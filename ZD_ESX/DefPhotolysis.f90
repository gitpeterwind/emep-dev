!> MODULE DefPhotolysis.f90 - A  crude version of the EMEP model's. Needs work.
!!****************************************************************************! 
!> Calculates Photolysis coefficients
!! Currently using "phux" method as documented in M. Kuhn et al., Atmos. 
!! Environ., 32, No.4, 693-709, 1998.
!! A second "MCM" method will be introduced soon. 
!! Compared to EMEP code, this routine calculates indices such as IDNO2, and
!! contains more descriptive structures.
!! 
!!------------------------------------------------------------------------------

     module DefPhotolysis
!-------------------------------------------------------------------------------

   implicit none
   private

 !/ subroutines: 

  public  :: setphotorates
  public  :: phuxj          !! As defined in Kuhn et al., AE, 1998
  public  :: InitPhotol     !! calls InitPhux or InitMCM (not yet)
  public  :: GetPhotol      !! 
  private :: InitPhux       !! 
  public  :: self_test      !! For testing this module

  private :: SetPhux        !! 
  private  :: ZenithAngle   !! For testing this module as stand-alone


 !/ parameters and variables: 

  character(len=10), private, parameter :: PhotolMethod = "Phux"
  real, private, parameter :: PI = 4.0*atan(1.0)

  integer, public, parameter :: &
       NRCPHOT      = 35  &! Number of pre-defined photolytic reactions
      ,MAXRCPHOT    = 41   ! Biggest index. TMP
     !EMEP   NRCPHOT      = 17   ! Number of pre-defined photolytic reactions

  integer, public, save :: nPhotol = 0
   
! In EMEP code, we use IDBO3 for O1D production and IDAO3 for OP production
!  Indices of photolysis rates as available from Phodis files:

!PHUX  !MCM coeffs: (not same as EMEP)
!PHUX    integer, public, parameter ::  &
!PHUX      IDAO3    =  1 , IDBO3    =  2 , IDNO2    =  4 , &
!PHUX!TMPTEST      IDAO3    =  2 , IDBO3    =  1 , IDNO2    =  4 , &
!PHUX      IDH2O2   =  3 , IDHNO3   =  8 , IDACH2O  =  11, &
!PHUX      IDBCH2O  = 12 , IDCH3CHO = 13 , IDCH3COX = 22 , &
!PHUX      !TMP IDCH3COY = -1 , IDHCOHCO = -1 , IDRCOHCO = 34 , &
!PHUX      IDCH3COY = -1 , IDRCOHCO = 34 , &
!PHUX      IDANO3   =  6 , IDN2O5   = -1 , IDCH3O2H = 41 , &
!PHUX      IDHO2NO2 = 16 , IDACETON = 17, &
!PHUX      IDAGLYOX = 31, IDBGLYOX = 32, IDCGLYOX = 33,&!RBMCM
!PHUX      IDBNO3   =  5                                !RBMCM, ca. 10% of IDNO3
!PHUX
!PHUX  !ESXTMP- FIX LATER
!PHUX    integer, public, parameter ::  &
!PHUX      IDNO3    =  IDANO3 !!!??????   IDANO3 -> NO2+ O3P
!PHUX
!PHUX    integer, public, parameter ::  IDRCOCHO  = IDRCOHCO ! Just tmp
!PHUX    integer, public, parameter ::  IDHCOHCO  = IDRCOHCO ! Just tmp

 !!--------------------------------------------------------------------------
 !! Phux system, from chemical comparison exercise of  Kuhn et al., AE, 1998
 !! which in turn was based upon Roeths program and RADM species.

  real, private, parameter :: MINYZ=-30.
  real, private, parameter :: EMINYZ=exp(MINYZ)
  integer, public, save :: IDO3_O3P, IDO3_O1D, &
     IDNO3, & ! SKIP: IDNO3_NO, IDNO3_NO2, &
     IDHONO, IDHCHO_H2, IDHCHO_H, IDCH3CHO, IDCH3O2H, IDH2O2, IDNO2, IDHNO3
  integer, public, save :: IDCH3COO2H,IDCH3COCH3,IDCHOCHO,IDRCOCHO,IDN2O5, IDMEK

  type, private :: Phx
    integer :: ind               ! index of e.g. IDNO2
    real :: x=0.0,y=0.0,z=0.0    ! x is photol rate for zenith angle 0
    character(len=50) :: reaction
  end type Phx

  type(Phx), public, dimension(NRCPHOT), save :: phux!=Phx()
 !!--------------------------------------------------------------------------


 !! From MCM

  type, private :: mcmdj
    integer :: ind
    real    :: L, M, N
    real    :: exj  ! example rate, mid latitudes, summer
  end type mcmdj

  type(mcmdj), dimension(35) :: dj = (/ &
   mcmdj(  1,   6.073E-05,   1.743,   0.474,   1.295E-06) &
  ,mcmdj(  2,   4.775E-04,   0.298,   0.080,   2.482E-04) &
  ,mcmdj(  3,   1.041E-05,   0.723,   0.279,   1.581E-06) &
  ,mcmdj(  4,   1.165E-02,   0.244,   0.267,   3.364E-03) &
  ,mcmdj(  5,   2.485E-02,   0.168,   0.108,   1.378E-02) &
  ,mcmdj(  6,   1.747E-01,   0.155,   0.125,   9.280E-02) &
  ,mcmdj(  7,   2.644E-03,   0.261,   0.288,   6.944E-04) &
  ,mcmdj(  8,   9.312E-07,   1.230,   0.307,   6.785E-08) &
  ,mcmdj( 11,   4.642E-05,   0.762,   0.353,   5.178E-06) &
  ,mcmdj( 12,   6.853E-05,   0.477,   0.323,   1.214E-05) &
  ,mcmdj( 13,   7.344E-06,   1.202,   0.417,   3.769E-07) &
  ,mcmdj( 14,   2.879E-05,   1.067,   0.358,   2.152E-06) &
  ,mcmdj( 15,   2.792E-05,   0.805,   0.338,   3.110E-06) &
  ,mcmdj( 16,   1.675E-05,   0.805,   0.338,   1.866E-06) &
  ,mcmdj( 17,   7.914E-05,   0.764,   0.364,   8.472E-06) &
  ,mcmdj( 18,   1.140E-05,   0.396,   0.298,   2.440E-06) &
  ,mcmdj( 19,   1.140E-05,   0.396,   0.298,   2.440E-06) &
  ,mcmdj( 21,   7.992E-07,   1.578,   0.271,   4.270E-08) &
  ,mcmdj( 22,   5.804E-06,   1.092,   0.377,   3.934E-07) &
  ,mcmdj( 23,   1.836E-05,   0.395,   0.296,   3.963E-06) &
  ,mcmdj( 24,   1.836E-05,   0.395,   0.296,   3.963E-06) &
  ,mcmdj( 31,   6.845E-05,   0.130,   0.201,   2.874E-05) &
  ,mcmdj( 32,   1.032E-05,   0.130,   0.201,   4.333E-06) &
  ,mcmdj( 33,   3.802E-05,   0.644,   0.312,   5.678E-06) &
  ,mcmdj( 34,   1.537E-04,   0.170,   0.208,   5.989E-05) &
  ,mcmdj( 35,   3.326E-04,   0.148,   0.215,   1.300E-04) &
  ,mcmdj( 41,   7.649E-06,   0.682,   0.279,   1.223E-06) &
  ,mcmdj( 51,   1.588E-06,   1.154,   0.318,   1.225E-07) &
  ,mcmdj( 52,   1.907E-06,   1.244,   0.335,   1.238E-07) &
  ,mcmdj( 53,   2.485E-06,   1.196,   0.328,   1.756E-07) &
  ,mcmdj( 54,   4.095E-06,   1.111,   0.316,   3.357E-07) &
  ,mcmdj( 55,   1.135E-05,   0.974,   0.309,   1.132E-06) &
  ,mcmdj( 56,   7.549E-06,   1.015,   0.324,   6.787E-07) &
  ,mcmdj( 57,   3.363E-06,   1.296,   0.322,   2.140E-07) &
  ,mcmdj( 61,   7.537E-04,   0.499,   0.266,   1.586E-04) &
  /)

 contains
 !>--------------------------------------------------------------------------

 !> Simple solar values for testing this module

    function ZenithAngle( doy, hr, lat ) result (theta)
      integer, intent(in) :: doy !! Day of year (1-366)
      real, intent(in) :: hr, lat
      real :: phi, lha, decl
      real :: theta

      lha = (1.0+hr/12)*PI         !! local hour angle
      phi = lat *PI/180.0    !! latitude in radians

     !! declination angle (varies daily)

      decl=-23.45*PI/180.0*cos(2.0*PI/365.0*(doy+10.0))

      theta = acos(cos(lha)*cos(decl)*cos(phi)+sin(decl)*sin(phi))

      !print *, "ZA ", theta*180./PI, decl*180./PI
    end function ZenithAngle

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine InitPhotol()
    call InitPhux()  !! Only option just now! Will code MCM later
  end subroutine InitPhotol
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function GetPhotol(idj,Za,debug_level) result (J)
    integer, intent(in) :: idj
    real,    intent(in) :: Za  ! Zenith angle, radians
    integer, intent(in) :: debug_level
    real :: J
    integer :: n, ind ! index in Phux array
      
    n=1
    do while( idj /= Phux(n)%ind )
      n=n+1
    end do
    if(debug_level>0) print *, "GETPHOTOL ", idj, n, &
       Phux(n)%ind, adjustl(Phux(n)%reaction)

    J= phuxj(n,Za)   ! Phux is only option just now

  end function  GetPhotol
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    function setphotorates( ind, cosZ, debug_level ) result(J) ! crude for ESX
      integer, intent(in) :: ind
      real, intent(in) :: cosZ !! cosine Solar zenith angle
      real :: J    !! Photorate
      integer, intent(in) :: debug_level
      real :: sun
      integer, dimension(MAXRCPHOT), save :: mapindex = -1
      integer :: imap

      MAPPING : if( mapindex(ind) < 1 ) then
        do imap = 1, size(dj%ind)
          if ( dj(imap)%ind == ind ) then
             mapindex(ind) = imap
             exit MAPPING
          end if
        end do
        print *, "DJ MAP FAILED", ind
        stop "DJ MAP FAILED"
      end if MAPPING

      imap = mapindex(ind)

      if( imap < 1 ) then
        print *, "DJ MAP NEG", ind
        stop "DJ MAP NEG"
      end if
        
      associate ( L=> dj(imap)%L, M=> dj(imap)%M, N => dj(imap)%N ) 

        if( cosZ > 1.0e-30 ) then
          J = L * cosZ**(M*exp(-n/cosZ))
        else
          J =  0.0
        end if
      end associate

    end function setphotorates
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine SetPhux(ind,x,y,z,reac)
    integer, intent(out) :: ind
    real, intent(in) :: x,y,z
    character(len=*), intent(in) :: reac
    nPhotol = nPhotol + 1
    phux(nPhotol) = phx( nPhotol, x,y,z,reac )
    !print *, "SetPhux ", nPhotol, trim(reac)
    ind = nPhotol
    
  end subroutine SetPhux
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> SUBROUTINE InitPhux
  !! Here we define X,Y,Z coefficients and reaction text. No need for all
  !! rates to be used.

  subroutine InitPhux()

    call SetPhux(IDNO2,     1.07e-2, 1.01319, 0.83330, "NO2->NO+O3P" ) 
    call SetPhux(IDO3_O1D,  3.22e-5, 4.45037, 0.83330, "O3->O1D") 
    call SetPhux(IDO3_O3P,  5.36e-4, 0.34764, 0.91030, "O3->O3P") 
    call SetPhux(IDHONO,    8.96e-4, 0.99438, 0.83925, "HONO->HO+NO") 
    call SetPhux(IDHNO3,    5.48e-7, 2.86922, 0.79561, "HNO3->HO+NO2") 
    !call SetPhux(IDNO3_NO,  2.74e-2, 0.26226, 0.92849, "NO3->NO+O2") 
    !call SetPhux(IDNO3_NO2, 2.73e-1, 0.29327, 0.92401, "NO3->NO2+O3P") 
    call SetPhux(IDNO3, 2.73e-1, 0.29327, 0.92401, "NO3->NO2+O3P")  ! from NO2 rate
    call SetPhux(IDH2O2,    7.78e-6, 1.91463, 0.79810, "H2O2->2OH" ) 
    call SetPhux(IDHCHO_H2, 4.92e-5, 1.60973, 0.80184, "HCHO->H2+CO") 
    call SetPhux(IDHCHO_H,  4.05e-5, 2.06917, 0.80267, "HCHO->H+HCO") !CO+2HO2
    call SetPhux(IDCH3CHO,  5.40e-6, 2.52915, 0.79722, "CH3CHO->CH3+CHO")
    call SetPhux(IDCH3O2H,  6.37e-6, 1.76570, 0.80004, "CH3OOH->CH3O+OH")  !CH3O2H sometimes
    call SetPhux(IDCH3COO2H,6.10e-9, 9.17009, 0.72585,  "CH3COO2H->products")
    call SetPhux(IDCH3COCH3,1.32e-5, 2.46350, 0.79768, "CH3COCH3->products")
    call SetPhux(IDCHOCHO,  3.11e-3, 0.55016, 0.88313, "CHOCHO->products") !GLYOX
    call SetPhux(IDRCOCHO,  1.85e-3, 0.57967, 0.87921, "CH3OCHO->products") !MGLYOX
    call SetPhux(IDN2O5,    3.79e-5, 1.70537, 0.80153, "N2O5->O3P+2NO2")
    call SetPhux(IDMEK,     1.36e-5, 2.59915, 0.80212, "CH3COC2H5->products")!CH3COX

  end subroutine InitPhux
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !! Calculate photolysis rate (1/s) for reaction idj and zenith angle Chi

  function phuxj(idj,Chi) result (J)

   integer, intent(in) :: idj
   real, intent(in) :: Chi  ! Zenith angle, radians
   real :: J
   real :: chiz, ychiz, eychiz
   real :: X,Y,Z

   X=phux(idj)%x
   Y=phux(idj)%y
   Z=phux(idj)%z

   chiz = Chi * Z
   if ( chiz < 1.57079632679489) then
     ychiz= Y*(1.0-(1.0/cos(ChiZ)))
     if ( ychiz > MINYZ ) then
       eychiz = exp(ychiz)
     else
       eychiz = EMINYZ
     endif
   else
       eychiz = EMINYZ
   end if
   J = X * eychiz
  end function phuxj

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !> Subroutine to run abovce calculations, and also demonstrate call-order

  subroutine self_test()
    real :: Za,cosZa  !! Solar zenith angle in radians, +cos(chi)
    real :: J, pi=4.0*atan(4.0)
    integer :: i

    Za =  ZenithAngle( doy=181, hr=12.0, lat=52.0 )
    cosZa = cos(Za)

    call InitPhotol()

    print *, "Zenith Angle ",Za*180./pi
    do i = 1, nPhotol
      print "(a,i3,1x,a20,g12.3)", "PHX ",i, &
       adjustl(phux(i)%reaction), phuxj(i , Za)
    end do
  end subroutine self_test

end module DefPhotolysis
!-----------------------------------------------------------------------------

!DSXprogram test_DefPhotolysis
!DSX  use DefPhotolysis
!DSX  implicit none
!DSX  call self_test()
!DSX!  print *, "NO2->NO+O3P"
!DSX!  print *, "MCM ", setphotorates( IDNO2, cosZa, 10 )
!DSX!!
!DSX!  print *, "NO3->NO+O2"
!DSX!  print *, "MCM ", setphotorates( IDBNO3, cosZa, 10 )
!DSX!
!DSX!  print *, "NO3->NO2+O3P"
!DSX!!  print *, "PHX ", Phux( 2.73e-1, 0.29327, 0.92401, Za)
!DSX!  print *, "MCM ", setphotorates( IDNO3, cosZa, 10 )
!DSX!
!DSX!  print *, "O3->O1D"
!DSX!  !print *, "MCMA", setphotorates( IDAO3, cosZa, 10 
!DSX!  print *, "MCMB", setphotorates( IDBO3, cosZa, 10 )
!DSX!
!DSX!  print *, "O3->O3P"
!DSXend program test_DefPhotolysis
