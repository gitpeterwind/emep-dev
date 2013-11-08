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
  use MicroMet_ml, only : AerRes !! TEST
  use ModelConstants, only : UNDEF_R

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
  real :: Ra_MO, za, Ka, n
  real, pointer, dimension(:) :: zlevs
  integer :: i, iz, nz
  integer :: nSL   !! index of surface layer ht. in zbnd
  logical :: first_call=.true.
  integer, parameter :: NFINE = 10
  real, dimension(NFINE) :: fineZ, fineKz  ! for fine-scale integration of Kz
  logical :: debug0

  debug0 = ( esx%debug_Zmet>0 .and. first_call ) 

  call CheckStop( L%Hmix < 0.0, "Hmix not set!" ) 
  
  nz = esx%nz
  zlevs => esx%zbnd(1:nz)

  L%hSL = 0.04 * L%Hmix          ! Height  of surface layer
  L%hSL = max( L%hSL, 2.0*Veg%h) ! Avoids matching problems in RSL
  nSL = count( zlevs < L%hSL+0.01 )

  if( debug0 ) write(*,*) "esx_Zmet HSL STUFF ", nz, &
       L%hSL, nSL, zlevs(nSL), size(zlevs)

  if( KzMethod == "constant" ) then
     if( first_call) print *, "KEYVALS", kwargs(1)
     Zmet%Kz = KeyValue( kwargs(1:1), "const" )

  else if( KzMethod == "power" ) then !test case
    za = KeyValue( kwargs, "za" )
    Ka = KeyValue( kwargs, "Ka" )
    n  = KeyValue( kwargs, "n" )
    Zmet%Kz =  def_kz_pow(zlevs, za, Ka, n )

  else if( KzMethod == "Leuning" ) then 

    Zmet(1:nz)%Kz = def_Kz_inc(zlevs,L%uStar,Veg%h,L%hSL,L%invL,Veg%dPerh,debug0 )

  end if

  !> Above surface layer we need to use other Kz values 

  if( KzMethod == "constant" ) then
    if( first_call) print *, "KzCONST", Zmet(3)%Kz
  else
    if ( L%invL > 0.0 ) then 

        Zmet(1:nz)%Kz2 =  JericevicKzE( zlevs,L%Hmix,L%ustar,0.0) 
        Zmet(nSL+1:nz)%Kz =  Zmet(nSL+1:nz)%Kz2

    else

        call O_BrienKz( L%hSL, L%Hmix, zlevs, L%ustar, L%invL, &
                    KzHs=Zmet(nSL)%Kz, KzHmix=0.0, Kz=Zmet(1:nz)%Kz3, &
                    debug_flag=debug0)
        Zmet(nSL+1:nz)%Kz =  Zmet(nSL+1:nz)%Kz3
    end if
  end if

  call print_Kz()

  ! Testing resistances. In general, Ra(z1,z2) = integral dz/Kz between z1 to z2
  ! within the surface layer (iz=1..nSL)
  ! MicroMet AerRes uses the analytical solution of standard M-O theory, so 
  ! knows nothing about the in-canopy (e.g. Leuning) type effects.

! if ( first_call ) then

    if(debug0) print "(a,2f8.2)", "AERO ustar,1/L ", L%ustar, L%invL
    do iz  = 1, nSL !-1

       fineZ(:) = (/ ( esx%z(iz)+(i-0.5)*esx%dzmid(iz)/NFINE, i=1,NFINE) /)

       fineKz(:) = def_Kz_inc(fineZ,L%uStar,Veg%h,L%hSL,L%invL,Veg%dPerh,debug0 )

       Ra_MO=AerRes( esx%z(iz), esx%z(iz+1), L%ustar, L%invL,0.4)

       Zmet(iz)%Ra = esx%dzmid(iz)*sum(1.0/fineKz)/NFINE

       if ( debug0 ) then
         print "(a,i3,f7.1,9(2x,a,2f10.2))", "AERO ", iz, esx%z(iz), &
           "Kz,dz: ", Zmet(iz)%Kz, esx%dzmid(iz) &
          ,"Ras: ", Ra_MO, esx%dzmid(iz)/Zmet(iz)%Kz &
          ,"from Int dz/Kz: ", Zmet(iz)%Ra
       end if

    end do

  first_call=.false.

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
     (/ "z    ", "dz   ", "dzmid","zbnd ", "zb/h ", "T(K) ", "Kz   ", "KzJG ", "KzOB " /) & 
        ,esx%z(1:nz)   & !coords
        ,reshape( (/ esx%dz(1:nz), esx%dzmid(1:nz), &
                     esx%zbnd(1:nz), esx%zbnd(1:nz)/Veg%h, &
             Zmet(1:nz)%tzK, &
             Zmet(1:nz)%Kz, Zmet(1:nz)%Kz2, Zmet(1:nz)%Kz3 /), (/ nz, 8 /) ),&
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
                                                               
