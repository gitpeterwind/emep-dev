
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

  use My_Derived_ml, only : f_2d, d_2d, D2_HMIX , &  !u7.4vg
        D2_FSTCF0, D2_FSTDF0, D2_FSTTC0, D2_FSTMC0, D2_FSTGR0, D2_FSTWH0,&
        D2_O3CF, D2_O3DF, D2_O3TC, D2_O3GR, D2_O3WH
  use Dates_ml,       only : date  
  use GenSpec_adv_ml  , only : NSPEC_ADV          &
                 ,IXADV_O3 &
                 ,IXADV_PAN , IXADV_NO , IXADV_NO2  &
                 ,IXADV_SO4,IXADV_SO2,IXADV_HNO3,IXADV_NH3                &
                 ,IXADV_aNH4, IXADV_aNO3,IXADV_pNO3    !6s                &
  use GenSpec_shl_ml, only:   & ! =>> IXSHL_xx
                IXSHL_OH,IXSHL_HO2
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
     NSITES_MAX =    51         & ! Max. no surface sites allowed
    ,FREQ_SITE  =    3          & ! Interval (hrs) between outputs
    ,NADV_SITE  =    3 &!NSPEC_ADV  & ! No. advected species (1 up to NSPEC_ADV)
    ,NSHL_SITE  =    1          & ! No. short-lived species
    ,NXTRA_SITE =    2            ! No. Misc. met. params  ( now T2)

   integer, public, parameter, dimension(NADV_SITE) :: &
    !SITE_ADV =  (/ (isite, isite=1,NADV_SITE) /)       ! Everything!!
    SITE_ADV =  (/ (isite, isite=1,3) /)       ! Everything!!

   integer, public, parameter, dimension(NSHL_SITE) :: &
    SITE_SHL =  (/ IXSHL_OH /)                          ! More limited!

  ! Extra parameters - need to be coded in Sites_ml also. So far
  ! we can choose from T2, or th (pot. temp.)

   character(len=10), public, parameter, dimension(NXTRA_SITE) :: &
     SITE_XTRA=  (/ "hmix    ", "th      "/) 
!    SITE_XTRA=  (/ "hmix    ", "Vg_ref  ", "Vg_1m   " /) 
!    SITE_XTRA=  (/ "th      ", "hmix    ", "Vg_ref  ", "Vg_1m   ", &
!   SITE_XTRA=  (/ "FSTCF0  ", "FSTDF0  ", "FSTTC0  ", "FSTMC0  ", &
!                  "FSTGR0  ", "FSTWH0  ", &  
!                  "O3CF    ", "O3DF    ", "O3TC    ", "O3GR    ", &
!                  "O3WH    " /) 

 !ds - rv1.6.12 - can access d_2d fields through index here, by
 !     setting "D2D" above and say D2_FSTCF0 here:

   integer,           public, parameter, dimension(NXTRA_SITE) :: &
    SITE_XTRA_INDEX=  (/     0,        0 /)     ! Height at mid-cell





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
     NSONDES_MAX =    51               &   ! Max. no sondes allowed
    ,FREQ_SONDE  =    12               &   ! Interval (hrs) between outputs
    ,NADV_SONDE  =    3                &   ! No.  advected species
    ,NSHL_SONDE  =    1                &   ! No. short-lived species
    ,NXTRA_SONDE =    1                    ! No. Misc. met. params  (now th)

   integer, public, parameter, dimension(NADV_SONDE) :: &
    SONDE_ADV =  (/ IXADV_O3, IXADV_NO, IXADV_NO2 /)
   integer, public, parameter, dimension(NSHL_SONDE) :: &
    SONDE_SHL =  (/ IXSHL_OH /)
   character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
    SONDE_XTRA=  (/ "z_mid" /)     ! Height at mid-cell

 !ds - rv1.6.12 - can access d_3d fields through index here, by
 !     setting "D3D" above and say D3_XKSIG12 here:

   integer,           public, parameter, dimension(NXTRA_SONDE) :: &
    SONDE_XTRA_INDEX=  (/     0 /)



   !==============================================================
   !/*** Hourly outputs   (from hourly_out routine) to print out
   !     concentrations every hour (or multiple: HOURLY_FREQ) for specified 
   !     sub-grid.
   !     Or even met. data (only temp2m specified so far  - others
   !     need change in hourly_out.f also).

   !----------------------------------------------------------------
   !ds rv1_8_2: Added possibility of multi-layer output. Specify
   ! NLEVELS_HOURLY here, and in hr_out defs use either:
   !
   !      ADVppbv to get surface concentrations (onyl relevant for 
   !              layer k=20 of course - gives meaningless  number f
   !               or higher levels.
   !Or,
   !      BCVppbv to get grid-centre concentrations (relevant for 
   !      all layers.
   !----------------------------------------------------------------

    integer, public, parameter :: NHOURLY_OUT =  2 ! No. outputs
    integer, public, parameter :: NLEVELS_HOURLY = 2 ! No. outputs
    integer, public, parameter :: FREQ_HOURLY = 3  ! 1 hours between outputs

    type, public:: Asc2D
         character(len=12):: name   ! Name (no spaces!)
         character(len=7) :: type   ! "ADVppbv" or "ADVugm3" or "SHLmcm3" 
         character(len=9) :: ofmt   ! Output format (e.g. es12.4)
         integer          :: spec   ! Species number in xn_adv or xn_shl array
                                    !ds u7.4vg .. or other arrays
         integer          :: ix1    ! bottom-left x
         integer          :: ix2    ! bottom-left y
         integer          :: iy1    ! upper-right x
         integer          :: iy2    ! upper-right y
         character(len=12) :: unit   ! Unit used 
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
        ISPEC_OUTBEG = 101  &
        ,JSPEC_OUTBEG = 51  &
        ,ISPEC_OUTEND = -100  &     ! Set negative here to exclude
        ,JSPEC_OUTEND = -50         ! Set negative here to exclude


   public :: set_output_defs

contains

 subroutine set_output_defs
   implicit none

   character(len=44) :: errmsg  ! Local error message
   integer           :: i       ! Loop index

   real                :: to_ug_S & ! conversion to ug of S
                         ,to_ug_N & ! conversion to ug of N
                         ,to_mgSIA& ! conversion to mg of N
                         ,to_ugSIA  ! conversion to ug of N 
   real, save :: m_s = 100.0 ! From cm/s to m/s

  ! introduce some integers to make specification of domain simpler
  ! and less error-prone. Numbers can be changed as desired.
   !merlininteger, save :: ix1 = 45, ix2 = 170, iy1=1, iy2 = 133 
   !merlin:
   !integer, save :: ix1 = 75, ix2 = 150, iy1=12, iy2 = 102 
   !integer, save :: ix1 = 110, ix2 = 115, iy1=60, iy2 = 60
   integer, save :: ix1 = 125, ix2 = 150, iy1=45, iy2 =  65  !! Athens/Thess. 


!jej - added 11/5/01 following Joffen's suggestion:

  to_ug_S = atwS*PPBINV/ATWAIR ! in output only accounting for Sulphur
  to_ug_N = atwN*PPBINV/ATWAIR ! in output only accounting for Nitrogen
  to_mgSIA= PPBINV/ATWAIR*1000.
  to_ugSIA= PPBINV/ATWAIR

 !/** Hourly outputs
 !    Note that the hourly output uses **lots** of disc space, so specify
 !    as few as you need and with as small format as possible (cf max value).

 ! ** REMEMBER : ADV species are mixing ratio !!
 ! ** REMEMBER : SHL species are in molecules/cm3, not mixing ratio !!
 ! ** REMEMBER : No spaces in name, except at end !!

  !**               name     type   
  !**                ofmt   ispec     ix1 ix2  iy1 iy2  unit conv    max

 !   include 'merlin_ixadv'
 hr_out(1)=  Asc2D("o3_50m", "BCVppbv", &
                  "(f9.4)",IXADV_o3, ix1,ix2,iy1,iy2, "ppbv",PPBINV,9000.0)
 hr_out(2)=  Asc2D("o3_3m", "ADVppbv", &
                  "(f9.4)",IXADV_o3, ix1,ix2,iy1,iy2, "ppbv",PPBINV,9000.0)

!  hr_out(1)=  Asc2D("Ozone", "ADVppbv", &
!                  "(f9.5)",IXADV_O3, ix1,ix2,iy1,iy2, "ppb",PPBINV,600.0)
!  hr_out(2)=  Asc2D("aNH4-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNH4, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(3)=  Asc2D("aNO3-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(4)=  Asc2D("SO4-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(5)=  Asc2D("pNO3-air", "ADVugm3", &
!                  "(f8.4)",IXADV_pNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,400.0)

!  hr_out(2)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_aNH4, 45, 170, 1, 133, "ug",to_ugSIA,600.0)
!  hr_out(3)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_aNO3, 45, 170, 1, 133, "ug",to_ugSIA,600.0)
!  hr_out(4)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_SO4, 45, 170, 1, 133, "ug",to_ugSIA,400.0)
!  hr_out(5)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_pNO3, 45, 170, 1, 133, "ug",to_ugSIA,400.0)
!
 ! Extra parameters - need to be coded in Sites_ml also. So far
 ! we can choose from T2, or th (pot. temp.)
 !ds u7.4vg - or from d_2d arrays.

  !**           type   ofmt   ispec    ix1 ix2  iy1 iy2  unit conv    max
  !hr_out(3)= &
  !   Asc2D("D2D", "(f6.1)",   D2_HMIX, ix1,ix2,iy1,iy2, "m",1.0   ,10000.0)

!ICP
!   hr_out(2)= Asc2D("Fst_TConif  ", "D2D", "(f8.5)",&
!                  D2_FSTCF0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(3)= Asc2D("Fst_TBroad  ", "D2D", "(f8.5)",&
!                  D2_FSTDF0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!IZ    hr_out(4)= Asc2D("Fst_MConif  ", "D2D", "(f8.5)",&
!IZ                   D2_FSTTC0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!IZ    hr_out(5)= Asc2D("Fst_MBroad  ", "D2D", "(f8.5)",&
!IZ                   D2_FSTMC0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(4)= Asc2D("Fst_Grass   ", "D2D", "(f8.5)",&
!                  D2_FSTGR0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(5)= Asc2D("Fst_Wheat   ", "D2D", "(f8.5)",&
!                  D2_FSTWH0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!IZ    hr_out(8)= Asc2D("O3__Wheat   ", "D2D", "(f7.3)",&
!IZ                   D2_O3WH,   ix1,ix2,iy1,iy2, "ppb       ", 1.0  ,900.0)
!IZ    hr_out(9)= Asc2D("O3__Beech   ", "D2D", "(f7.3)",&
!IZ                   D2_O3DF,   ix1,ix2,iy1,iy2, "ppb       ", 1.0  ,900.0)
!IZ    hr_out(10)= Asc2D("O3__Conif   ", "D2D", "(f7.3)",&
!IZ                   D2_O3CF,   ix1,ix2,iy1,iy2, "ppb       ", 1.0  ,900.0)

 !/** theta is in deg.K
 !hr_out(1)=  Asc2D("T2_C",   "T2_C   ", &
 !                "(f5.1)",     -99, ix1,ix2,iy1,iy2, "degC",1.0   ,100.0)
 !hr_out(2)=  Asc2D("Precip", "PRECIP ", &
 !                "(f11.7)",    -99, ix1,ix2,iy1,iy2, "mm/hr",1.0,  200.0)
 !hr_out(3)=  Asc2D("Idir",   "Idirect", &
 !                "(f5.1)",    -99, ix1,ix2,iy1,iy2, "umole/m2/s",1.0, 1000.0)
 !hr_out(4)=  Asc2D("Idif",   "Idiffus", &
 !                "(f5.1)",    -99, ix1,ix2,iy1,iy2, "umole/m2/s",1.0, 1000.0)





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

     wanted_dates_bi(1) = date(-1,1,1,0,0)
     wanted_dates_bi(2) = date(-1,1,1,3,0)
     wanted_dates_bi(3) = date(-1,1,1,6,0)

 end subroutine set_output_defs
end module My_Outputs_ml




