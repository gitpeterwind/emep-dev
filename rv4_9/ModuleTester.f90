program tester
  use AeroFunctions,          only: testAeroFunctions     =>self_test
! use AeroFunctions,          only: testAeroFunctions     =>self_test_fracs
  use GridAllocate_ml,        only: testGridAllocate      =>self_test
  use KeyValueTypes,          only: testKeyValueTypes     =>self_test
  use TimeDate_ExtraUtil_ml,  only: testTimeDate_ExtraUtil=>self_test
  use Io_Progs_ml,            only: testIo_Progs          =>self_test
! use DefPhotolysis,          only: testDefPhotolysis     =>self_test_fracs
  use SmallUtils_ml,          only: testSmallUtils        =>self_test
  implicit none
  integer :: i
  character(len=64) :: ml
          
  do i = 1, iargc()
    call getarg(i, ml)
    write(*,*) "Testing: ",trim(ml)
    select case(ml)
    case("AeroFunctions","AeroFunctions.f90")
      call testAeroFunctions()
    case("GridAllocate_ml","GridAllocate_ml.f90")
      call testGridAllocate()
    case("KeyValueTypes","KeyValueTypes.f90")
      call testKeyValueTypes()
    case("TimeDate_ExtraUtil_ml","TimeDate_ExtraUtil_ml.f90")
      call testTimeDate_ExtraUtil()
    case("Io_Progs_ml","Io_Progs_ml.f90")
      call testIo_Progs()
!   case("DefPhotolysis","DefPhotolysis.f90")
!     call testDefPhotolysis()
    case("SmallUtils_ml","SmallUtils_ml.f90")
      call testSmallUtils()
    endselect
  enddo
endprogram tester
