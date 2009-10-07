module OwnDataTypes_ml
  implicit none

    public :: print_dep_type
    integer, public, parameter :: TXTLEN_DERIV = 24

  ! Contains some user-defined data-types, and routine associated
  ! with these. Collecting them here will
  ! avoid some dependencies, and shorten some My type modules.
  !
  ! Deriv from old Derived_ml
  ! Dep_type from old My_Derived_ml
  ! (should be able to merge these two in next step)
  ! VBST from SOA_ml

   !================== 
    type, public:: Deriv
       character(len=10) :: class ! Type of data, e.g. ADV or VOC
       logical  :: avg      ! True => average data (divide by nav at end),
                            !     else accumulate over run period
       integer  :: index    ! index in concentation array, or other
       real     :: scale    ! Scaling factor
       logical  :: rho      ! True when scale is ug (N or S)
       logical  :: inst     ! True when instantaneous values needed
       logical  :: year     ! True when yearly averages wanted
       logical  :: month    ! True when monthly averages wanted
       logical  :: day      ! True when daily averages wanted
       character(len=TXTLEN_DERIV) :: name ! Name of the variable (for netCDF output)
       character(len=10) :: unit ! Unit (writen in netCDF output)
    end type Deriv

 ! Store indices of ecossytem-species depositions
 ! dep_type( name, LC, index, f2d, class, label, txt, scale, atw, units )
 !            x     d      d    d   a10    a10   a10     f    i    a10
  type, public :: Dep_type
     character(len=TXTLEN_DERIV) :: name ! e.g. DDEP_SO2_m2Conif
     integer :: LC            ! Index of Receiver land-cover (one 
                              ! of Dep_Receivers)
     integer :: Index         ! e.g. index in e.g. xn_adv arrays, or Y for AFstY
     integer :: f2d           ! index in f_2d arrays
     character(len=10) :: class !  "Mosaic"
     character(len=10) :: label !  e.g. "VG", "Rns"
     character(len=10) :: txt ! text where needed, e.g. "Conif"
     real    :: scale         !  e.g. use 100.0 to get cm/s
     integer :: atw           ! atomic weight where needed
     character(len=10) :: units ! e.g.  mgN/m2
  end type 

  !+ Defines SOA, NONVOL and VBS params 

  type, public :: VBST
       integer     :: index    ! just for clarity
       real        :: CiStar   ! ug/m3
      !real        :: Tref    ! Assumed 300
       real        :: DeltaH   ! kJ/mole
  end type VBST



contains
 !=========================================================================
subroutine print_dep_type(w)
  type(dep_type), intent(in) :: w  ! wanted

  write(6,*) "Prints dep_type ========================="
  write(6,"(a,a)")      "Name   :", w%name
  write(6,"(a,i3)")     "LC     :", w%LC
  write(6,"(a,i3)")     "index  :", w%index
  write(6,"(a,i3)")     "f2d    :", w%f2d
  write(6,"(a,a10)")    "class  :", w%class
  write(6,"(a,a10)")    "label  :", w%label
  write(6,"(a,a10)")    "txt    :", w%txt
  write(6,"(a,es10.3)") "scale  :", w%scale
  write(6,"(a,i3)")     "atw    :", w%atw
  write(6,"(a,a10)")    "units  :", w%units
end subroutine print_dep_type

 !=========================================================================
end module OwnDataTypes_ml
