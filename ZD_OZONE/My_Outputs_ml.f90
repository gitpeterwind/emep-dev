
module  My_Outputs_ml

  !/ Allows user to specify which species are output to various
  !  ascii and binary output files.
  !
  ! From Steffen's Output_defs_ml, re-coded by ds to apply to site,
  ! sonde and hourly ouytputs, and to use "type" arrays.
  !
  ! Sites  - surface sites,     to sites.out
  ! Sondes - vertical profiles, to sondes.out
  ! Hourly - ascii output of selected species, selcted domain
  ! Restri - Full 3-D output of all species, selected domain

  use Dates_ml,       only : date  
  use GenSpec_adv_ml, only : NSPEC_ADV          &
        ,IXADV_O3 , IXADV_PAN , IXADV_CO , IXADV_NO , IXADV_NO2  &
        ,IXADV_SO4,IXADV_SO2,IXADV_HNO3,IXADV_NH3                &
        ,IXADV_HCHO,IXADV_CH3CHO, IXADV_NO3                      &
        ,IXADV_N2O5 , IXADV_AMSU, IXADV_AMNI    !6s                &
!6s         ,IXADV_SOA, IXADV_ASOA, IXADV_BSOA
  use GenSpec_shl_ml, only:   & ! =>> IXSHL_xx
                IXSHL_OH,IXSHL_HO2,  NSPEC_SHL
  !u1 use GenSpec_maps_ml, only:    MAP_ADV2TOT
  use GenChemicals_ml , only: species
  use ModelConstants_ml, only: PPBINV, PPTINV, ATWAIR, atwS, atwN
  use Par_ml,            only: me, NPROC   ! for gc_abort
  implicit none
  private


   !/*** Site outputs   (used in Sites_ml)
   !==============================================================
   !     Specify the species to be output to the sites.out file
   ! For met params we have no simple index, so we use characters.
   ! These must be defined in Sites_ml.f90.

   integer, private :: isite              ! To assign arrays, if needed
   integer, public, parameter :: &
     NSITES_MAX =    30         & ! Max. no surface sites allowed
    ,FREQ_SITE  =    1          & ! Interval (hrs) between outputs
    ,NADV_SITE  =    NSPEC_ADV  & ! No. advected species (1 up to NSPEC_ADV)
    ,NSHL_SITE  =    NSPEC_SHL  & ! No. short-lived species
    ,NXTRA_SITE =    1            ! No. Misc. met. params  ( now T2)

   integer, public, parameter, dimension(NADV_SITE) :: &
    !jej SITE_ADV =  (/ IXADV_NO, IXADV_NO2, IXADV_O3 /) 
    SITE_ADV =  (/ (isite, isite=1,NADV_SITE) /)       ! Everything!!

   integer, public, parameter, dimension(NSHL_SITE) :: &
!!    SITE_SHL =  (/ IXSHL_OH /)                          ! More limited!
    SITE_SHL =  (/ (isite, isite=1,NSHL_SITE) /)       ! TROTREP everything!!

  ! Extra parameters - need to be coded in Sites_ml also. So far
  ! we can choose from T2, or th (pot. temp.)

   character(len=10), public, parameter, dimension(NXTRA_SITE) :: &
    SITE_XTRA=  (/ "th" /)



   !/*** Aircraft outputs   (used in Polinat_ml)
   !==============================================================
   !     Specify the species to be output by Polinat for aircraft flight tracks

   integer, public, parameter :: &
     NFLIGHT_MAX =    10               &   ! Max. no sondes allowed
    ,FREQ_FLIGHT =    12               &   ! Interval (hrs) between outputs
    ,NADV_FLIGHT =    1                    ! No.  advected species

   integer, public, parameter, dimension(NADV_FLIGHT) :: &
    FLIGHT_ADV =  (/ IXADV_O3 /)



   !/*** Sonde outputs   (used in Sites_ml)
   !==============================================================
   !     Specify the species to be output to the sondes.out file
   !  We typically deal with fewer species for sonde output than
   !  surface sites, so we use a different method to specify.
   ! For met params we have no simple index, so we use characters.
   ! These must be defined in Sites_ml.f90.

   integer, public, parameter :: &
     NSONDES_MAX =    30               &   ! Max. no sondes allowed
    ,FREQ_SONDE  =    12               &   ! Interval (hrs) between outputs
    ,NADV_SONDE  =    4                &   ! No.  advected species
    ,NSHL_SONDE  =    1                &   ! No. short-lived species
    ,NXTRA_SONDE =    1                    ! No. Misc. met. params  (now th)

   integer, public, parameter, dimension(NADV_SONDE) :: &
    SONDE_ADV =  (/ IXADV_O3, IXADV_NO, IXADV_NO2, IXADV_CO /)
   integer, public, parameter, dimension(NSHL_SONDE) :: &
    SONDE_SHL =  (/ IXSHL_OH /)
   character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
    SONDE_XTRA=  (/ "z_mid" /)     ! Height at mid-cell



   !==============================================================
   !/*** Hourly outputs   (from hourly_out routine) to print out
   !     concentrations every hour (or multiple: HOURLY_FREQ) for specified 
   !     sub-grid.
   !     Or even met. data (only temp2m specified so far  - others
   !     need change in hourly_out.f also).

    integer, public, parameter :: NHOURLY_OUT = 1  ! No. outputs
    integer, public, parameter :: FREQ_HOURLY = 1  ! 3 hours between outputs

    type, public:: Asc2D
         character(len=3) :: type   ! "ADV" or "SHL" 
        character(len=12) :: ofmt   ! Output format (e.g. es12.4)
         integer          :: spec   ! Species number in xn_adv or xn_shl array
         integer          :: ix1    ! bottom-left x
         integer          :: ix2    ! bottom-left y
         integer          :: iy1    ! upper-right x
         integer          :: iy2    ! upper-right y
         character(len=4) :: unit   ! Unit used 
         real             :: unitconv   !  conv. factor
         real             :: max    ! Max allowed value for output
    end type Asc2D
     
    type(Asc2D), public, dimension(NHOURLY_OUT) :: hr_out  ! Set below


  !/** wanted binary dates... specify days for which full binary
  !    output is wanted. Replaces the hard-coding which was
  !    in wrtchem:

     integer, public, parameter :: NBDATES = 3  
     type(date), public, save, dimension(NBDATES) :: wanted_dates_bi 



!/*** Restricted outputs   (calls restri_out outine) 
!==============================================================
!  Output of all species, all levels, for restricted domain.
!  parameters for the domain, which should be written - specify!!!!
!    if you do not want to print on restricted domain:
!    set upper bound (END) smaller than lower bound (BEG)

!  integer, public, parameter ::  &
!       ISPEC_OUTBEG = ISMBEG  &
!       ,JSPEC_OUTBEG = JSMBEG  &
!       ,ISPEC_OUTEND = ISMBEG+GIMAX-1  &
!       ,JSPEC_OUTEND = JSMBEG+GJMAX-1
  integer, public, parameter ::  &
        ISPEC_OUTBEG = 70  &
        ,JSPEC_OUTBEG = 25  &
        ,ISPEC_OUTEND = -125  &     ! Set negative here to exclude
        ,JSPEC_OUTEND = -80         ! Set negative here to exclude

!   integer, public, parameter :: NPARUPP = NSPEC_ADV ! 39

! Defines integers which represent output species (former mod_def):

!6c  integer, public, parameter ::  & 
!6c       NSPEC_ADD_ADV = 1        !& ! additional ps advection, if = 1
!       , NPARADV = NSPEC_ADV + NSPEC_ADD_ADV  ! total species to be advected


   public :: set_output_defs

contains

 subroutine set_output_defs
   use GenSpec_adv_ml, only:  IXADV_O3  ! for Hourly outputs
   use GenSpec_shl_ml, only:  IXSHL_OH   ! for Hourly outputs
   implicit none

   character(len=30) :: errmsg  ! Local error message
   integer           :: i       ! Loop index

   real                :: to_ug_S & ! conversion to ug of S
                         ,to_ug_N   ! conversion to ug of N
 
 
!jej - added 11/5/01 following Joffen's suggestion:

  to_ug_S = atwS*PPBINV/ATWAIR ! in output only accounting for Sulphur
  to_ug_N = atwN*PPBINV/ATWAIR ! in output only accounting for Nitrogen


 !/** Hourly outputs
 !    Note that the hourly output uses **lots** of disc space, so specify
 !    as few as you need and with as small format as possible (cf max value).

 ! ** REMEMBER : ADV species are mixing ratio !!
 ! ** REMEMBER : SHL species are in molecules/cm3, not mixing ratio !!

  !**           type   ofmt   ispec    ix1 ix2  iy1 iy2  unit conv    max

  hr_out(1)= &
    Asc2D("ADV", "(f9.5)",IXADV_O3, 55, 150, 10, 100, "ppb",PPBINV,600.0)
 ! hr_out(2)= &
 !   Asc2D("ADV", "(f8.4)",IXADV_NO, 55, 150, 10, 100, "ppb",PPBINV,600.0)
 !! hr_out(3)= &
 !   Asc2D("ADV", "(f8.4)",IXADV_NO2, 55, 150, 10, 100, "ppb",PPBINV,600.0)
 ! hr_out(4)= &
 !   Asc2D("ADV", "(f8.2)",IXADV_CO, 55, 150, 10, 100, "ppb",PPBINV,4000.0)
!    Asc2D("ADV", "(g9.3)",IXADV_SOA, 70, 150, 10, 100, "ppb",PPBINV,400.0)
!  hr_out(3)= &
!    Asc2D("ADV", "(g9.3)",IXADV_ASOA, 70, 150, 10, 100, "ppb",PPBINV,400.0)
!  hr_out(4)= &
!    Asc2D("ADV", "(g9.3)",IXADV_BSOA, 70, 150, 10, 100, "ppb",PPBINV,400.0)
!!!hr_out(2)= &
!!  Asc2D("SHL", "(f8.4)",IXSHL_OH, 87, 150, 51, 75, "10^5",1.0e-5,999.0)

 ! Extra parameters - need to be coded in Sites_ml also. So far
 ! we can choose from T2, or th (pot. temp.)

 !/** theta is in deg.K
!  hr_out(5)= &
!     Asc2D("th ", "(f5.1)",     -99, 70, 150, 10, 100, "degK",1.0   ,400.0)

 !/**Asc2D("T2 ", "(f5.1)",     -99, 87, 110, 51, 75, "degC",1.0   ,100.0)


  !/** Consistency checks

   do i = 1, NHOURLY_OUT

       ! We use ix1 to see if the array has been set.

        if ( hr_out(i)%ix1 < 1 .or.  hr_out(i)%ix1 > 999 ) then
            errmsg = "Failed consistency check in set_output_defs"
            print *,  errmsg, "Hourly: ", i, "Nhourly: ", NHOURLY_OUT
            call gc_abort(me,NPROC,errmsg)
        end if
   end do

  !/** wanted binary dates... specify days for which full binary
  !    output is wanted. Replaces the hard-coding which was
  !    in wrtchem:

     wanted_dates_bi(1) = date(-1,6,6,12,0)
     wanted_dates_bi(2) = date(-1,6,8,12,0)
     wanted_dates_bi(3) = date(-1,6,14,12,0)

 end subroutine set_output_defs
end module My_Outputs_ml
