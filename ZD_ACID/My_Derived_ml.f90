
!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module contains the definitions and function used to get
  ! "derived" fields, such as accumulated precipitation or sulphate, 
  ! daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as binary fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  ! 
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined 
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml. 
  !  
  ! These 2 modules combine Steffen's rewriting of the output routines
  ! (IOU notation), elements of chemint_mach, Bud_ml, etc.
  ! Re-coded to use only 2 types of data (d_2 and d_3)
  ! and F90 types for Deriv by ds, Sept. 2001. 

  !
  ! ds, 15/10/01
  !---------------------------------------------------------------------------
 
use GenSpec_adv_ml        ! Use IXADV_ indices...
use GenSpec_tot_ml,  only : SO4 !, HCHO, CH3CHO     !  For mol. wts.
use GenChemicals_ml, only : species               !  For mol. wts.
use ModelConstants_ml, only : atwS, atwN, ATWAIR  &
                        , KMAX_MID &  ! =>  z dimension
                        , PPBINV  &   !   1.0e9
                        , MFAC    &   ! converts roa (kg/m3 to M, molec/cm3)
                        , current_date

!6c - needed for derived fields
use Chemfields_ml, only : xn_adv, xn_shl, cfac
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use Met_ml,        only : z_bnd, roa    ! 6c REM: zeta
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area 
use PhysicalConstants_ml,  only : PI
use Radiation_ml,  only :  zen
implicit none
private

 public  :: Set_My_Derived 
 public  :: My_DerivFunc 

 private :: acc_sulphate         ! Sums sulphate column
! private :: aot_calc             ! Calculates daylight AOTs
 private :: misc_xn              ! Miscelleaneous Sums and fractions of xn_adv


  ! Define first the 4 possible former output types
  ! corresponding to instantaneous,year,month,day

   integer, public, parameter ::  & 
        IOU_INST=1, IOU_YEAR=2, IOU_MON=3, IOU_DAY=4

  ! 6c: Replace old NIOUTP with separate 2d and 3d dimensions:
  ! The 2-d and 3-d fields use the above as a time-dimension. We define
  ! LENOUTxD according to how fine resolution we want on output. For 2d
  ! fields we use daily outputs. For the big 3d fields, monthly output
  ! is sufficient.

   integer, public, parameter ::  LENOUT2D = 4  ! Allows INST..DAY for 2d fields
   integer, public, parameter ::  LENOUT3D = 3  ! Allows INST..MON for 3d fields

  ! ***  ds 26/9/2001 - compressed arrays of steffen into new types
  !       for Deriv3D and Deriv2D. Values now set in subroutine
  !       Set_My_Derived

    type, public:: Deriv    ! Could be private ??
       integer  :: code     ! Identifier for DNMI/xfelt (was ID6OUT_DERIV_3D)
       character(len=7) :: class ! Type of data, e.g. ADV or VOC
       logical  :: avg      ! True => average data (divide by nav at end), 
                            !     else accumulate over run period
       integer  :: index    ! index in concentation array, or other
       real     :: scale    ! Scaling factor          ! (was SCALOUT_DERIV_3D)
       logical  :: rho      ! True when scale is ug (N or S)
       logical  :: inst     ! True when instantaneous values needed
       logical  :: year     ! True when yearly averages wanted
       logical  :: month    ! True when monthly averages wanted
       logical  :: day      ! True when daily averages wanted
    end type Deriv

   logical, private, parameter :: T = .true., F = .false. ! shorthands only

   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.

   integer, public, parameter ::  & 
        NWDEP = 4       &   ! Number of 2D deposition fields
       ,WDEP_PREC  = 1  &   ! sum rainfall, was IPRDEP
       ,WDEP_SOX   = 2  &   ! sum of sulphur        (was xwdep (IXWD_SOX)
       ,WDEP_OXN   = 3  &   ! sum oxidised nitrogen (was  xddep(IXWD_HNO3))
       ,WDEP_RDN   = 4      ! sum reduced  nitrogen (was  xddep(IXWD_NH3))

   integer, public, parameter ::  & 
        NDDEP = 3       &   ! Number of 2D deposition fields
       ,DDEP_SOX   = 1  &   ! sum of sulphur        (was xwdep (IXDD_SOX)
       ,DDEP_OXN   = 2  &   ! sum oxidised nitrogen (was  xddep(IXDD_HNO3))
       ,DDEP_RDN   = 3      ! sum reduced  nitrogen (was  xddep(IXDD_NH3))

   integer, public, parameter ::  & 
        NDERIV_2D = 15 &   ! Number of 2D derived fields
       ,D2_ACCSU  = 1  &   ! was NUM_ACCSU 
       ,D2_SO2    = 2  &   ! was xnsurf(so2)
       ,D2_SO4    = 3  &   ! was xnsurf(..)
       ,D2_HNO3   = 4  &   ! was xnsurf(..)
       ,D2_PAN    = 5  &   ! was xnsurf(..)
       ,D2_NH3    = 6  &   ! was xnsurf(..)
       ,D2_NO2    = 7  &   ! was xnsurf(..)
       ,D2_AMSU   =8   &   ! was xnsurf(..)
!hf hmix
       ,D2_HMIX   = 9     &
       ,D2_HMIX00 = 10    &!mixing height at 00
       ,D2_HMIX12 = 11     !mixing height at 12
   integer, public, parameter ::  & 
        D2_SOX      = 12 &  ! Sum of sulphates, 
       ,D2_OXN      = 13  &  ! Total nitrates (HNO3 + part. NITRATE) 
       ,D2_REDN     = 14  &  ! Annonia + ammonium 
       ,D2_FRNIT    = 15     ! (part. nitrate)/(tot. nitrate)

   integer, public, parameter ::  & 
        NDERIV_3D = 6    & ! Number of 3D derived fields
       ,D3_SO2    = 1   &  ! was xnav(no2) array
       ,D3_NO2    = 2   &
       ,D3_NO    = 3    &
       ,D3_PAN = 4      &
       ,D3_XKSIG00= 5   &
       ,D3_XKSIG12=6
   ! We put definitions in f_2d, f_3d, and  data into d_2d, d_3d:

    type(Deriv), public, dimension(NWDEP),     save :: f_wdep!wet dep
    type(Deriv), public, dimension(NDDEP),     save :: f_ddep!dry dep
    type(Deriv), public, dimension(NDERIV_2D), save :: f_2d  !other deriv
    type(Deriv), public, dimension(NDERIV_3D), save :: f_3d  !other deriv

   ! Note - previous versions did not have the LENOUT2D dimension
   ! for wet and dry deposition. Why not?  Are annual or daily
   ! depositions never printed? Since I prefer to keep all 2d
   ! fields as similar as posisble, I have kept this dimension
   ! for now - ds

    real, save,  public :: &
      wdep( NWDEP    ,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !wet dep
      ddep( NDDEP    ,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !dry dep
      d_2d( NDERIV_2D,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !other deriv
      d_3d( NDERIV_3D,MAXLIMAX, MAXLJMAX, KMAX_MID, LENOUT3D )

    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine Set_My_Derived()

   ! Set the parameters for the derived parameters, including the codes
   ! used by DNMI/xfelt and scaling factors. (The scaling factors may
   ! be changed later in Derived_ml.
   
    real          :: sf           ! scaling factor
    real, save    :: ugS = atwS*PPBINV/ATWAIR
    real, save    :: ugN = atwN*PPBINV/ATWAIR


   !-- Deposition fields

   ! Deriv type has fields:     code class   avg? ind scale rho  Inst Yr  Mn   Day   
    f_wdep( WDEP_PREC ) = Deriv( 561, "PREC ", F, -1, 1.0,   F  , F  ,  T , T ,  T )
    f_wdep( WDEP_SOX  ) = Deriv( 541, "WDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )
    f_wdep( WDEP_OXN  ) = Deriv( 542, "WDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )
    f_wdep( WDEP_RDN  ) = Deriv( 543, "WDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )

    f_ddep( DDEP_SOX  ) = Deriv( 521, "DDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )
    f_ddep( DDEP_OXN  ) = Deriv( 522, "DDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )
    f_ddep( DDEP_RDN  ) = Deriv( 523, "DDEP ", F, -1, 1.0e6, F  , F  ,  T , T ,  T )

   !-- 2-D fields - the complex ones

   ! Deriv type has fields:  code class  avg? ind scale rho  Inst  Yr  Mn   Day
!    f_2d( D2_AOT40) = Deriv( 608, "AOT  ", F, 40, 1.0,   F  , F  ,  T , T ,  F  )
!    f_2d( D2_AOT60) = Deriv( 609, "AOT  ", F, 60, 1.0,   F  , F  ,  T , T ,  F  )

     ! To give same units as MACHO....  (multiplied with roa in layers?? ==> rho "false" )
    sf = species ( SO4 )%molwt * PPBINV /ATWAIR
    f_2d( D2_ACCSU) = Deriv( 611, "ACCSU", T, -1, sf ,   F  , F  ,  T , T ,  F  )

   ! -- simple advected species

   ! Deriv type has fields:     code class   avg? ind scale rho  Inst Yr  Mn   Day   
    f_2d( D2_SO2 ) = Deriv( 601, "ADV  ", T, IXADV_SO2, ugS, T  , F ,  T , T ,  T )
    f_2d( D2_SO4 ) = Deriv( 620, "ADV  ", T, IXADV_SO4, ugS, T  , F ,  T , T ,  T )
    f_2d( D2_HNO3) = Deriv( 621, "ADV  ", T, IXADV_HNO3,ugN, T  , F ,  T , T ,  T )
    f_2d( D2_PAN ) = Deriv( 604, "ADV  ", T, IXADV_PAN, ugN, T  , F ,  T , T ,  T )
    f_2d( D2_NH3 ) = Deriv( 623, "ADV  ", T, IXADV_NH3, ugN, T  , F ,  T , T ,  T )
    f_2d( D2_NO2 ) = Deriv( 606, "ADV  ", T, IXADV_NO2, ugN, T  , F ,  T , T ,  T )
    f_2d( D2_AMSU) = Deriv( 619, "ADV  ", T,IXADV_AMSU, ugS, T  , F ,  T , T ,  T )
!hf hmix
    f_2d( D2_HMIX)   = Deriv( 468, "HMIX    ", T      ,0.0 , 1.0, T  , F ,  T , T ,  T )
    f_2d( D2_HMIX00) = Deriv( 469, "HMIX00  ", T      ,0.0 , 1.0, T  , F ,  T , T ,  T )
    f_2d( D2_HMIX12) = Deriv( 470, "HMIX12  ", T      ,0.0 , 1.0, T  , F ,  T , T ,  T )
!    f_2d( D2_O3  ) = Deriv( 607, "ADV  ", T, IXADV_O3 ,PPBINV, F  , F , T , T , T )
!    f_2d( D2_CO  ) = Deriv( 612, "ADV  ", T, IXADV_CO ,PPBINV, F  , F , T , T , T )

   ! --  time-averages - here 8-16 , as used in MACHO

!    sf = species ( HCHO )%molwt * PPBINV /ATWAIR
!    f_2d( D2T_HCHO)  =Deriv( 613,"TADV ", T, IXADV_HCHO  ,sf , T , F , T , T , T)

!    sf = species ( CH3CHO )%molwt * PPBINV /ATWAIR
!    f_2d( D2T_CH3CHO)=Deriv( 614,"TADV ", T, IXADV_CH3CHO,sf , T , F , T , T , T)

!    f_2d( D2T_VOC  ) =Deriv( 610,"TVOC ", T,   -1    ,PPBINV , F , F , T , T , T)


   ! -- miscellaneous user-defined functions

   ! Deriv type has fields:     code class   avg? ind scale rho  Inst Yr  Mn   Day  
    f_2d( D2_SOX   ) =    Deriv( 602,"TSO4 ", T,   -1  ,ugS , T , F , T , T , T)

    f_2d( D2_OXN   ) =    Deriv( 603,"TOXN ", T,   -1  ,ugN , T , F , T , T , T)

    f_2d( D2_REDN   ) =   Deriv( 605,"TRDN ", T,   -1  ,ugN , T , F , T , T , T)

    f_2d( D2_FRNIT   ) =  Deriv( 624,"FRNIT", T,   -1  ,1.0 , F , F , T , T , T)

   !-- 3-D fields

    f_3d( D3_SO2  ) = Deriv( 401, "ADV  ", T, IXADV_SO2 , ugS , T , T , T , T , F )
    f_3d( D3_NO2 ) = Deriv( 406, "ADV  ", T, IXADV_NO2, PPBINV , F , T , T , T , F )
    f_3d( D3_NO ) = Deriv( 407, "ADV  ", T,  IXADV_NO , PPBINV , F , T , T , T , F )
    f_3d( D3_PAN ) = Deriv( 408, "ADV  ", T,  IXADV_PAN , PPBINV , F , T , T , T , F )
    f_3d( D3_XKSIG00 ) = Deriv( 409, "XKSIG00  ", T,  0. , 1. , F , T , T , T , F )
    f_3d( D3_XKSIG12 ) = Deriv( 410, "XKSIG12  ", T,  0. , 1. , F , T , T , T , F )
  end subroutine Set_My_Derived
 !=========================================================================
  subroutine My_DerivFunc( n, class , timefrac, density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml. For example, we need the 
    ! index for IXADV_O3 for AOTs, and this might not be available in the model
    ! we are running (a PM2.5 model for example), so it is better to define 
    ! this function here.

  integer, intent(in) :: n           ! number of output data field
  character(len=*), intent(in)    :: class       ! Class of data
  real, intent(in)    :: timefrac    ! Timestep as frationof hour, dt/3600

  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density     
! density = 1 ( or = roa when unit ug)

  select case ( class )

      case ( "ACCSU" )

          call acc_sulphate(n)

!      case ( "AOT" )

!           call aot_calc( n, timefrac )

      case ( "TSO4", "TOXN", "TRDN", "FRNIT"  )

!!print *, "Calling misc_xn for ", class
           call misc_xn( n, class, density )

      case  default

            print *, "WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
     end select
  


  end subroutine My_DerivFunc
 !=========================================================================

  subroutine acc_sulphate( n )

    !/--  adds up sulphate column 

    integer, intent(in) :: n               ! Index for output field

    forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST) +  &
             ( xn_adv(IXADV_SO4,i,j,KMAX_MID) +    &
               xn_adv(IXADV_AMSU,i,j,KMAX_MID) ) * &
               (z_bnd(i,j,k) - z_bnd(i,j,k+1)) *   &
                 roa(i,j,k,1)*1.0e9
    end forall
  end subroutine acc_sulphate
 !=========================================================================
!hf aot out
 !=========================================================================

  subroutine misc_xn( n , class, density)
    integer, intent(in) :: n       ! Index for output field
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  
! density = 1 ( or = roa when unit ug)


    !/--  adds up sulphate, nitrate, or whatever is defined

    select case ( class )

    case ( "TSO4" ) 
      forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = &
               ( xn_adv(IXADV_SO4,i,j,KMAX_MID) * cfac(IXADV_SO4,i,j)  &
               + xn_adv(IXADV_AMSU,i,j,KMAX_MID) * cfac(IXADV_AMSU,i,j) ) &
               * density(i,j)
      end forall


    case ( "TOXN" )
      forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_AMNI,i,j,KMAX_MID) * cfac(IXADV_AMNI,i,j) )&
!hf Changed for MADE version
!              + xn_adv(IXADV_NITRATE,i,j,KMAX_MID) * cfac(IXADV_NITRATE,i,j)) &
              * density(i,j)
      end forall


    case ( "TRDN" )
      forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = &
               ( xn_adv(IXADV_NH3,i,j,KMAX_MID) * cfac(IXADV_NH3,i,j)    &
              +  xn_adv(IXADV_AMNI,i,j,KMAX_MID) * cfac(IXADV_AMNI,i,j)  &
          + 1.5* xn_adv(IXADV_AMSU,i,j,KMAX_MID) * cfac(IXADV_AMSU,i,j)) & 
               * density(i,j)
      end forall

!hf What is this????
    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = &
                xn_adv(IXADV_AMNI,i,j,KMAX_MID) * cfac(IXADV_AMNI,i,j)  &
            / (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j)   &
            +  xn_adv(IXADV_AMNI,i,j,KMAX_MID) * cfac(IXADV_AMNI,i,j) )    
      end forall
!!d_2d( n, i,j,IOU_INST) +  &

!!print *, "misc_xn  After  :" , n, d_2d(n,2,2,IOU_INST)
    end select
!    print *, "total conc. ", xn_adv(IXADV_SO4,3,3,KMAX_MID),    &
!                 xn_adv(IXADV_AMSU,3,3,KMAX_MID) 
  end subroutine misc_xn
 !=========================================================================
end module My_Derived_ml


