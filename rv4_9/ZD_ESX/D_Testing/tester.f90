! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!UNCOMMENT TO TEST JUST THIS MODULE:
program tester
  use esx_Zgrid_ml,  only : test_Zgrid, nzlev ! TESTING only
  use esx_Zveg_ml,   only : test_Zveg ! TESTING only
  integer :: ionum = 6

  ! If file wanted, use other ionum for TestingLAI.txt output
  !open(ionum,file="TestingLAI.txt")

  call  test_Zgrid(ionum) 
  print *, "Into Zveg nzlev=", nzlev
  call  test_Zveg(ionum) 

  !close(ionum)
end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
