
!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module specifies the "derived" fields, such as accumulated precipitation 
  ! or sulphate, daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as netCDF fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  ! 
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined 
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml. 
  !  
  ! ds, modified 16/12/2003:
  !   Major changes. Only text strings used here to define wanted data
  !   All data field characteristics should be defined in Derived_ml, e.g.
  !   in f_2d arrays. 
  !   Derived fields such as d_2d only exist in Derived_ml, so are
  !   accessed here through subroutine calls - using just the (i,j) part
  !   of the bigger d_2d arrays
  !---------------------------------------------------------------------------
 
use GenSpec_adv_ml        ! Use IXADV_ indices...
use GenSpec_shl_ml        ! Use IXSHL_ indices...
use GenSpec_tot_ml,  only : SO4, HCHO, CH3CHO  &   !  For mol. wts.
                           ,aNO3, pNO3, aNH4, PM25, PMCO &
                           ,SSfi, SSco  !SeaS
use GenChemicals_ml, only : species               !  For mol. wts.
use ModelConstants_ml, only : atwS, atwN, ATWAIR  &
                        , KMAX_MID &  ! =>  z dimension
                        , PPBINV  &   !   1.0e9
                        , MFAC    &   ! converts roa (kg/m3 to M, molec/cm3)
                        , current_date

use Chemfields_ml, only : xn_adv, xn_shl, cfac
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use Met_ml,        only : z_bnd, roa    ! 6c REM: zeta
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area 
!ds use PhysicalConstants_ml,  only : PI
!ds use Radiation_ml,  only :  zen
implicit none
private

 !ds public  :: Set_My_Derived 
 public  :: My_DerivFunc 

 !BUGGY? private :: acc_sulphate         ! Sums sulphate column
 !ds moved private :: aot_calc             ! Calculates daylight AOTs
 private :: misc_xn   &          ! Miscelleaneous Sums and fractions of xn_adv
           ,pm_calc              ! Miscelleaneous PM's


   character(len=8),  public ,parameter :: model='ZD_OZONE'

   integer, public, parameter :: NOUTPUT_ABS_HEIGHTS = 6

   real, public, parameter, dimension(NOUTPUT_ABS_HEIGHTS) :: &
          OUTPUT_ABS_HEIGHTS = &
             (/ 0.0, 1.0, 3.0, 5.0, 10.0, 20.0 /)  ! Above ground
            ! Note - values below d+z0m will be set to d+z0m

   logical, private, parameter :: T = .true., F = .false. ! shorthands only

   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.
   !  Factor 1.0e6 converts from kg/m2/a to mg/m2/a

  !ds 24/3/2004:
  !dsNSR - new system for distinguishing source-receptor (SR) stuff from model
  !        evaluation.  The SR runs should use as few as possible outputs
  !        to keep CPU and disc-requirements down. We define first then the
  !        minimum list of outputs for use in SR, then define an extra list
  !        of parameters needed in model evaluation, or even for the base-case
  !        of SR runs.



  !============ parameters for source-receptor modelling: ===================!

  !ds 24/3/2004:
  ! This logical variable is not used here, but in Derived_ml, we set all IOU_DAY
  ! false if SOURCE_RECPTOR = .true.. We don't (generally) want daily outputs
  ! for SR runs.

    logical, public, parameter :: SOURCE_RECEPTOR = .false.


  ! Then we have some standard SR outputs.. ( a bit longer than necessary right now)

    integer, public, parameter :: NSR_2D = 37

    character(len=12), public, parameter, dimension(NSR_2D) :: &
  D2_SR = (/ &
!
!    Particles: components
       "D2_SO4      ","D2_aNO3     ","D2_pNO3     ","D2_aNH4     " &
      ,"D2_PPM25    ","D2_PPMco    ","D2_PM25_H2O " &
!
!    Particles: sums
      ,"D2_SIA      ","D2_PM25     ","D2_PM10     ","D2_PMco     " &
      ,"D2_SS       ","D2_tNO3     " &
!
!    Ozone and AOTs
      ,"D2_O3       ","D2_MAXO3    " &
      ,"D2_AOT30    ","D2_AOT40    ","D2_AOT60    " &
      ,"D2_EUAOT30WH","D2_EUAOT30DF","D2_EUAOT40WH","D2_EUAOT40DF" &
      ,"D2_UNAOT30WH","D2_UNAOT30DF","D2_UNAOT40WH","D2_UNAOT40DF" &
!
!    NOy-type sums 
      ,"D2_NO2      ","D2_OXN      ","D2_NOX      ","D2_NOZ      " &
!
!    Ecosystem - fluxes:
      ,"D2_FSTDF00  ","D2_FSTDF16  ","D2_FSTWH00  ","D2_FSTWH30  " &
      ,"D2_FSTWH60  ","D2_O3DF     ","D2_O3WH     " &
  /)

  !============ Extra parameters for model evaluation: ===================!

    integer, private, parameter :: NEXTRA_2D =  7
    character(len=12), public, parameter, dimension(NEXTRA_2D) :: &
  D2_EXTRA = (/ &
       "D2_SO2      ","D2_HNO3     ","D2_NH3      ","D2_NO       "&
      ,"D2_REDN     ","D2_SSfi     ","D2_SSco     " &
  /)


   !============ Choose here (Comment out EXTRA stuff for SR runs ============!
   ! Number of other 2D derived fields used:

   integer, public, parameter :: NDERIV_2D = &
           !NSR_2D                   ! SOURCE_RECEPTOR = .true.
           NSR_2D + NEXTRA_2D       ! SOURCE_RECEPTOR = .false.

   character(len=12), public, parameter, dimension(NDERIV_2D) :: &
      D2_USED = (/  D2_SR &
                   ,D2_EXTRA /)   ! SOURCE_RECEPTOR = .false.
                   !         /)   ! SOURCE_RECEPTOR = .true.


!----------------------
! Less often needed:
 !exc  "D2_CO     ","D2T_HCHO  ","D2T_CH3CHO","D2_VOC    ",
 !exc ,"D2_O3CF   ","D2_O3TC   ","D2_O3GR   ","D2_ACCSU  ",
 !"D2_FRNIT  ","D2_MAXOH  ","D2_HMIX   ","D2_HMIX00 ","D2_HMIX12 " &
 !exc  "D2_PAN    ","D2_AOT20    " /)

   !=== NEW MY_DERIVED SYSTEM ======================================
   !ds 16/12/2003

   ! Define number of fields for each type:

  integer, public, parameter :: &
        NWDEP     =  4   &  !
       ,NDDEP     = 15   &  ! wetlands and water now removed
!dsNSR ,NDERIV_2D = 45   &  ! Number of other 2D derived fields used  !water&!SeaS
       ,NDERIV_3D =  0      ! Number of 3D derived fields

  ! then use character arrays to specify which are used.

   character(len=9), public, parameter, dimension(NWDEP) :: &
       WDEP_USED = (/ "WDEP_PREC", "WDEP_SOX ", "WDEP_OXN ", &
                      "WDEP_RDN " /)   ! WDEP_PM not used

  !ds waters and wetlands removed:

   character(len=10), public, parameter, dimension(NDDEP) :: &
     DDEP_USED = (/  &
        "DDEP_SOX  ","DDEP_OXN  ","DDEP_RDN  "  &
       ,"DDEP_OXSCF","DDEP_OXSDF","DDEP_OXSCR","DDEP_OXSSN"  &
       ,"DDEP_OXNCF","DDEP_OXNDF","DDEP_OXNCR","DDEP_OXNSN"  &
       ,"DDEP_RDNCF","DDEP_RDNDF","DDEP_RDNCR","DDEP_RDNSN"  &
     /)    !  "DDEP_PM   " not needed?


     character(len=13), public, parameter, dimension(0:NDERIV_3D) :: &
!     character(len=13), public, parameter, dimension(NDERIV_3D) :: &
       D3_USED = (/ "DUMMY" /) ! Dummy value if empty
!ds    D3_USED = (/"D3_O3        ","D3_OH        ", "D3_CH3COO2   ","D3_PHNO3     "&
!ds               ,"D3_H2O2      "   &
!ds               ,"D3_MAXOH     ", "D3_MAXCH3COO2" /) ! Dummy value if empty

    !ds - lines defining ddep, wdep, etc. moved to Derived_ml
!====

    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine My_DerivFunc( e_2d, n, class , timefrac, density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml. For example, we need the 
    ! index for IXADV_O3 for AOTs, and this might not be available in the model
    ! we are running (a PM2.5 model for example), so it is better to define 
    ! this function here.

  real, dimension(:,:), intent(inout) :: e_2d  !  (i,j) 2-d extract of d_2d
  integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
  character(len=*), intent(in)    :: class       ! Class of data
  real, intent(in)    :: timefrac    ! Timestep as frationof hour, dt/3600

  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density     
! density = 1 ( or = roa when unit ug)

  select case ( class )

      !BUGGY? case ( "ACCSU" )

      !BUGGY?     call acc_sulphate(n)

      !ds MOVED to Derived_ml
      !MOVED case ( "AOT" )  ! ds added AFSTy
      !MOVED      call aot_calc( e_2d, n, ndef, timefrac )

      !ds rv1_9_17 case ( "TSO4", "TOXN", "TRDN", "FRNIT", "tNO3 "    )
      !ds rv1_9_32 case ( "TOXN", "TRDN", "FRNIT", "tNO3 " , "SSalt"   )
      case ( "NOX", "NOZ", "TOXN", "TRDN", "FRNIT", "tNO3 " , "SSalt"   )

!!print *, "Calling misc_xn for ", class
           call misc_xn( e_2d, n, class, density )

      case ( "SIA", "PM10", "PM25", "PMco" )

          call pm_calc(e_2d, n, class,  density)

      case  default

            print *, "WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
     end select
  


  end subroutine My_DerivFunc
 !=========================================================================
!BUGGY?
!BUGGY?  subroutine acc_sulphate( acc_2d )
!BUGGY?
!BUGGY?   !/--  adds up sulphate column 
!BUGGY?
!BUGGY?   !ds BUGGY - started to change this, then noticed that "k" wasn't set
!BUGGY?   !           and KMAX_MID used in xn_adv. Just comment out for now.
!BUGGY?   !ds NEW***
!BUGGY?   !ds - now pass in i,j part of d_2d array with 
!BUGGY?   !     say call acc_sulphate(d_2d(n,:,:,IOU_INST)
!BUGGY?
!BUGGY?   real, dimension(:,:), intent(inout) :: acc_2d   ! Extract of 2d field
!BUGGY?
!BUGGY?    forall ( i=1:limax, j=1:ljmax )
!BUGGY?        acc_2d( i,j) = acc_2d( i,j) +  &
!BUGGY?              xn_adv(IXADV_SO4,i,j,KMAX_MID)*     &
!BUGGY?              (z_bnd(i,j,k) - z_bnd(i,j,k+1)) *   &
!BUGGY?                 roa(i,j,k,1)*1.0e9
!BUGGY?    end forall
!BUGGY?  end subroutine acc_sulphate
 !=========================================================================
   !MOVED: ds
   !subroutine aot_calc MOVED to Derived_ml
 !=========================================================================
 !=========================================================================

  subroutine pm_calc( pm_2d, n, class, density )

    !/--  calulates PM10 = SIA + PPM2.5 + PPMco

    real, dimension(:,:), intent(inout) :: pm_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  

    select case ( class )

    case ( "SIA" ) 

      forall ( i=1:limax, j=1:ljmax )
        !ds d_2d( n, i,j,IOU_INST) = &
        pm_2d( i,j) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j))& 
         * density(i,j)
      end forall

    case ( "PM25" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt *cfac(IXADV_SSfi,i,j))&  !SeaS
         * density(i,j)
      end forall

    case ( "PMco" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSco,i,j,KMAX_MID) *species(SSco)%molwt *cfac(IXADV_SSco,i,j))&  !SeaS
          * density(i,j)
      end forall

    case ( "PM10" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt*cfac(IXADV_SO4,i,j)   &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt*cfac(IXADV_SSfi,i,j) & !SeaS
         + xn_adv(IXADV_SSco,i,j,KMAX_MID)*species(SSco)%molwt*cfac(IXADV_SSco,i,j))& !SeaS
         * density(i,j)
      end forall

    end select

  end subroutine  pm_calc
 !=========================================================================

!=========================================================================

  subroutine misc_xn( e_2d, n, class, density)
    real, dimension(:,:), intent(inout) :: e_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  
! density = 1 ( or = roa when unit ug)


    !/--  adds up sulphate, nitrate, or whatever is defined

    select case ( class )

    !case ( "TSO4" )  !ds rv1_9_17 - relic of old system - remove!
    !  forall ( i=1:limax, j=1:ljmax )
    !      !ds d_2d( n, i,j,IOU_INST) = &
    !      e_2d( i,j ) = xn_adv(IXADV_SO4,i,j,KMAX_MID) * cfac(IXADV_SO4,i,j)  &
    !                   * density(i,j)
    !  end forall


    case ( "TOXN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
              * density(i,j)
      end forall

!ds 31/3/04 .. added new groupngs for SR: NOX, NOZ
! Shoudl allow calculation of NOy = NOx + NOz

    case ( "NOX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_NO,i,j,KMAX_MID) &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) &
              ) * density(i,j)
      end forall

    case ( "NOZ" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) &
              + xn_adv(IXADV_PAN,i,j,KMAX_MID) * cfac(IXADV_PAN,i,j) &
              + xn_adv(IXADV_MPAN,i,j,KMAX_MID) * cfac(IXADV_MPAN,i,j) &
              + xn_adv(IXADV_NO3,i,j,KMAX_MID) &
              + 2.0* xn_adv(IXADV_N2O5,i,j,KMAX_MID) &
              + xn_adv(IXADV_ISNI,i,j,KMAX_MID) &
              ) * density(i,j)
      end forall


    case ( "TRDN" )
      forall ( i=1:limax, j=1:ljmax )
          !ds d_2d( n, i,j,IOU_INST) = &
          e_2d( i,j ) = &
               ( xn_adv(IXADV_NH3,i,j,KMAX_MID) * cfac(IXADV_NH3,i,j)    &
              +  xn_adv(IXADV_aNH4,i,j,KMAX_MID) * cfac(IXADV_aNH4,i,j))  &
               * density(i,j)
      end forall


    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
             ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)  &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
            /max(1E-80, (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j))   &
            +  xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)    &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j))
      end forall

    case ( "tNO3" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d(  i,j ) = &
              ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) )&
              * density(i,j)
      end forall

!SeaS
    case ( "SSalt" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_SSfi,i,j,KMAX_MID) * cfac(IXADV_SSfi,i,j) &
              + xn_adv(IXADV_SSco,i,j,KMAX_MID) * cfac(IXADV_SSco,i,j) )&
              * density(i,j)
      end forall

    end select
  end subroutine misc_xn
 !=========================================================================
end module My_Derived_ml


