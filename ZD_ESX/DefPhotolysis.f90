!> MODULE DefPhotolysis.f90 - A  crude version of the EMEP model's. Needs work.
!!****************************************************************************! 

!> Calculates Photolysis coefficients
!!------------------------------------------------------------------------------

     module DefPhotolysis
!-------------------------------------------------------------------------------

   implicit none
   private

   integer, public, parameter :: &
       NRCPHOT      = 35  &! Number of pre-defined photolytic reactions
      ,MAXRCPHOT    = 41   ! Biggest index. TMP
     !EMEP   NRCPHOT      = 17   ! Number of pre-defined photolytic reactions
   
!  Indices of photolysis rates as available from Phodis files:

  !MCM coeffs: (not same as EMEP)
    integer, public, parameter ::  &
      IDAO3    =  2 , IDBO3    =  1 , IDNO2    =  4 , &
      IDH2O2   =  3 , IDHNO3   =  8 , IDACH2O  =  11, &
      IDBCH2O  = 12 , IDCH3CHO = 13 , IDCH3COX = 22 , &
      !TMP IDCH3COY = -1 , IDHCOHCO = -1 , IDRCOHCO = 34 , &
      IDCH3COY = -1 , IDRCOHCO = 34 , &
      IDANO3   =  6 , IDN2O5   = -1 , IDCH3O2H = 41 , &
      IDHO2NO2 = 16 , IDACETON = 17, &
      IDAGLYOX = 31, IDBGLYOX = 32, IDCGLYOX = 33,&!RBMCM
      IDBNO3   =  5                                !RBMCM, ca. 10% of IDNO3

  !ESXTMP- FIX LATER
    integer, public, parameter ::  &
      IDNO3    =  IDBNO3 !!!??????
    real, private, save :: pi = 4.0*atan(1.0)

    integer, public, parameter ::  IDRCOCHO  = IDRCOHCO ! Just tmp
    integer, public, parameter ::  IDHCOHCO  = IDRCOHCO ! Just tmp

 !/ subroutines: 

  public :: setphotorates

  type, private :: mcmdj
    integer :: ind
    real    :: L
    real    :: M
    real    :: N
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
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine setphotorates( ind, time, J,  debug_level ) ! crude for ESX
      integer, intent(in) :: ind
      real, intent(in) :: time
      real, intent(out) :: J    !! Photorate
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
        

      !just to get started, use sun only to scale mcm rates

      sun=cos( (time/3600.0 - 12)/24.0 * 2 * pi )

      if( sun > 0.0 ) then
          J = sun * dj(imap)%exj
      else
          J =  0.0
      end if

      !> Simple debug=1 is treated external to this routine. Only
      !! activate for debug>=2
      if( debug_level >= 2 .and. ind == IDNO2 ) &
          write(*,"(a,2i3,f8.3,es12.3)") "RCPHOT-NO2", ind, imap, sun, J

    end subroutine setphotorates
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 end module DefPhotolysis
