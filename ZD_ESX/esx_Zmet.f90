!>  MODULE esx_Zmet
!!
!!  Definition of ESX-specific meteorology for z-layers
!!  We assign some initial values to simplify testing 
!   (Can/should be reset by real programmes!)

module esx_Zmet
  use CheckStops,    only: CheckStop
  use Io_Routines,   only: writedata
  use KeyValueTypes, only: KeyValue, KeyValReal
  use esx_Variables, only: L=> Loc, esx, Zmet
  use esx_Zveg, only: Veg
  use Kz_ml                  !! e.g. def_Kz, JericevicKz , O_BrienKz

  implicit none
  private

  public ::  set_esxZmet  ! allocates arrays, sets initial(fake) values
  public ::  def_kz
  public ::  print_Kz
  public ::  test_Zmet  ! 

  logical, public, save :: first_call = .true.


 contains

  !----------------------------------------------------------------------------
  !> init_Zmet - mainly for testing. Initialises some met from local surface
  !! data. VERY CRUDE !!! Should have some lapse rates and other d/dz changes

  subroutine set_esxZmet()
    integer :: nz
    nz = esx%nz

   !> Initialise for test

    if ( Zmet(1)%tzK < 1 )   then ! Assume nothing set.
      if( esx%debug_Zmet>0) print "(a,i4, 2f7.2,es10.3)", "INIT_ZMET from Local, t,rh,p =>", &
        nz, L%t2,L%rh,L%psurf
      Zmet(1:nz)%tzK = L%t2   
      Zmet(1:nz)%rh  = L%rh
      Zmet(1:nz)%Pa  = L%psurf
    end if

    call def_Kz( esx%Kz_method, esx%Kz_kwargs )

  end subroutine set_esxZmet

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> def_kz calculates Kz values for given z-levels.
!! A variety of methods are available, with parameters provided through
!! the key-word argument methodology
 
subroutine def_Kz(KzMethod, kwargs)

  character(len=*), intent(in) :: KzMethod
  type(KeyValReal), dimension(:) :: kwargs
  real :: za, Ka, n, kappa = 0.4   ! CHANGE kappa later
  real, pointer, dimension(:) :: zlevs
  integer :: nz, nSL   !! index of surface layer ht. in zbnd


  call CheckStop( L%Hmix < 0.0, "Hmix not set!" ) 
  
  nz = esx%nz
  zlevs => esx%zbnd(1:nz)

  L%hSL = 0.04 * L%Hmix ! Height  of surface layer
  nSL = count( zlevs < L%hSL+0.01 )

  if( esx%debug_Zmet>0) write(*,*) "esx_Zmet HSL STUFF ", nz, &
       L%hSL, nSL, zlevs(nSL), size(zlevs)

  if( KzMethod == "constant" ) then
     print *, "KEYVALS", kwargs(1)
     Zmet%Kz = KeyValue( kwargs(1:1), "const" )

  else if( KzMethod == "power" ) then !test case
    za = KeyValue( kwargs, "za" )
    Ka = KeyValue( kwargs, "Ka" )
    n  = KeyValue( kwargs, "n" )
    Zmet%Kz =  def_kz_pow(zlevs, za, Ka, n )

  else if( KzMethod == "Leuning" ) then 

    za = KeyValue( kwargs, "za" )
    Ka = KeyValue( kwargs, "Ka" )
    n  = KeyValue( kwargs, "n" )

    Zmet(1:nz)%Kz = def_Kz_inc(zlevs,L%uStar,Veg%hVeg,L%hSL,L%invL,kappa,Veg%dPerh)

  end if

  !> Above surface layer we need to use other Kz values 

  if ( L%invL > 0.0 ) then 

    Zmet(1:nz)%Kz2 =  JericevicKzE( zlevs,L%Hmix,L%ustar,0.0) 
    Zmet(nSL+1:nz)%Kz =  Zmet(nSL+1:nz)%Kz2

  else

    call O_BrienKz( L%hSL, L%Hmix, zlevs, L%ustar, L%invL, &
                    KzHs=Zmet(nSL)%Kz, KzHmix=0.0, Kz=Zmet(1:nz)%Kz3, &
                    debug_flag=esx%debug_Zmet>0 )
    Zmet(nSL+1:nz)%Kz =  Zmet(nSL+1:nz)%Kz3
  end if
  

  call print_Kz()

end subroutine def_kz
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!> print_Kz provides neat table of Kz values.
!! Can be called with ionum=6 to get standard output, otherwise LogKz.txt.
!! Assumed to only be called once, otherwise LogKz.txt overwritten (?)

 subroutine print_Kz()
   character(len=100) :: header
   integer :: nz
   header = "Kz method = "// trim(esx%Kz_method)

   ! need to reshape 1-D arrays into 2-D for writedata:

   nz = esx%nz
   call writedata("LogKz", &
     (/ "z    ", "dz   ", "dzmid","zbnd ", "T(K) ", "Kz   ", "KzJG ", "KzOB " /) & 
        ,esx%z(1:nz)   & !coords
        ,reshape( (/ esx%dz(1:nz), esx%dzmid(1:nz), &
                     esx%zbnd(1:nz), Zmet(1:nz)%tzK, &
             Zmet(1:nz)%Kz, Zmet(1:nz)%Kz2, Zmet(1:nz)%Kz3 /), (/ nz, 7 /) ),&
      header )

 end subroutine  print_Kz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !----------------------------------------------------------------------------
  !> Simple test. Not much should go wrong though (but be careful with z as
  !! mid-point or grid boundary)

  subroutine test_zmet(ionum)
    integer, intent(in) :: ionum    ! Unit number for output file
    real, dimension(5) :: ztest = (/1.0, 4.0, 8.0, 16.0, 100.0 /)
    integer :: iz
    esx%nz = 5

    call set_esxZmet()

    write(ionum,"(a)") "Test Met nz=", size(ztest)
    do iz = 1, size(ztest)
      write(ionum,"(a,i3,f6.1,4f12.3)") "TestMet", iz,ztest(iz),Zmet(iz)%tzK, Zmet(iz)%rh,Zmet(iz)%Pa
    end do
    call  def_Kz( esx%Kz_method, esx%Kz_kwargs)
    call  print_Kz() 
  end subroutine test_Zmet
  !----------------------------------------------------------------------------
end module esx_Zmet

!NB: TESTING PROG AT END IF WANTED
!SKIP  !> config_Zmet 
!SKIP
!SKIP  subroutine config_Zmet(io,writelog)
!SKIP    integer, intent(in) :: io
!SKIP    logical, intent(in) :: writelog
!SKIP    integer :: ilog
!SKIP
!SKIP    namelist /esxMeteo_config/Zmet
!SKIP    print *, "INTO ZMET===================================="
!SKIP    rewind(io)
!SKIP    read (io, nml=esxMeteo_config)
!SKIP    if(  writelog ) then
!SKIP      open(newunit=ilog,file="LogConfig.Zmet")
!SKIP      print *, "CONFIG ZMET====================================", ilog
!SKIP      write(ilog,"(a)") "CONFIG ZMET===================================="
!SKIP      write(ilog,nml=esxMeteo_config)
!SKIP      write(*,nml=esxMeteo_config)
!SKIP      close(ilog)
!SKIP    end if
!SKIP  end subroutine config_Zmet

!UNCOMMENT TO TEST JUST THIS MODULE:
!DSX program tester
!DSX   use LocalVariables, only : L
!DSX   use esx_Zgrid, only : init_Zgrid, test_Zgrid
!DSX   use esx_Zmet,  only : test_Zmet, config_Zmet, set_esxZmet
!DSX   integer :: ionum
!DSX    
!DSX    ! Some initial values of ambient temperature
!DSX    L%t2 = 298.0; L%psurf=1.0e5; L%rh=0.6
!DSX    open(newunit=ionum,file="config_esx.nml")
!DSX    call config_Zmet(ionum,writelog=.false.)
!DSX    close(ionum)
!DSX  print *, "DONE CONFIG"
!DSX
!   If file wanted, use other ionum for TestingMet.txt output
!
!DSX   open(newunit=ionum,file="TestingMet.txt")
!DSX    call  test_Zgrid(ionum) 
!DSX    call   set_esxZmet() 
!DSX    call  test_Zmet(ionum) 
!DSX 
!DSX   close(ionum)
!DSX end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
                                                               
