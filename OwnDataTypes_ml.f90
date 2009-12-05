module OwnDataTypes_ml
  implicit none

    public :: print_deriv_type
    integer, public, parameter :: TXTLEN_DERIV = 24
    integer, public, parameter :: TXTLEN_SHORT  = 12

  ! Contains some user-defined data-types, and routine associated
  ! with these. Collecting them here will
  ! avoid some dependencies, and shorten some My type modules.
  !
  ! Deriv used i My_Derived and  Derived_ml
  ! VBST from SOA_ml

   !================== 

    type, public:: Deriv
       character(len=TXTLEN_DERIV) :: name ! e.g. DDEP_SO2_m2Conif
       character(len=TXTLEN_SHORT) :: class ! Type of data, e.g. ADV or Mosaic
       character(len=TXTLEN_SHORT) :: subclass !  e.g. "VG", "Rns"
       character(len=TXTLEN_SHORT) :: txt ! text where needed, e.g. "Conif"
       character(len=TXTLEN_SHORT) :: unit ! writen in netCDF output
       integer  :: index    ! index in concentation array, or other
       integer :: f2d           ! index in f_2d arrays
       integer :: LC            ! Index of Receiver land-cover (one 
       real    :: XYCL          ! Threshold or CL, e.f. AOTx or AFstY
       real    :: scale         !  e.g. use 100.0 to get cm/s
       real :: dt_scale     ! used only if we need a factor on dt_advec,
                            ! otherwise set to zero.
       logical  :: avg      ! True => average data (divide by nav at end),
                            !     else accumulate over run period
       logical  :: rho      ! True when scale is ug (N or S)
       logical  :: inst     ! True when instantaneous values needed
       logical  :: year     ! True when yearly averages wanted
       logical  :: month    ! True when monthly averages wanted
       logical  :: day      ! True when daily averages wanted
                              ! of Dep_Receivers)
       integer :: atw           ! atomic weight where needed
  end type 

  !================== 
  !+ Defines SOA, NONVOL and VBS params 

  type, public :: VBST
       integer     :: index    ! just for clarity
       real        :: CiStar   ! ug/m3
      !real        :: Tref    ! Assumed 300
       real        :: DeltaH   ! kJ/mole
  end type VBST

contains
 !=========================================================================

subroutine print_Deriv_type(w)
  type(Deriv), intent(in) :: w  ! wanted

  write(6,*) "Prints Deriv type ========================="
  write(6,"(a,a)")      "Name   :", w%name
  write(6,"(a,a10)")    "class  :", w%class
  write(6,"(a,a10)")    "subclass  :", w%subclass
  write(6,"(a,a10)")    "txt    :", w%txt   
  write(6,"(a,a10)")    "units  :", w%unit
  write(6,"(a,i3)")     "LC     :", w%LC
  write(6,"(a,i3)")     "index  :", w%index
  write(6,"(a,i3)")     "f2d    :", w%f2d
  write(6,"(a,a10)")    "txt    :", w%txt
  write(6,"(a,es10.3)") "scale  :", w%scale
  write(6,"(a,i3)")     "atw    :", w%atw
end subroutine print_Deriv_type

 !=========================================================================
end module OwnDataTypes_ml
