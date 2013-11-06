module esx_ChemRun
  use ChemSpecs
  use ChemRates
  use ChemSolver
  use DefPhotolysis   !! Need IDNO2etc as well as setphotorates, NRCPHOT
  use esx_Variables, only: esx, Zmet, Loc
  use PhysicalConstants_ml, only : DEG2RAD
  use ZchemData, only : xChem, Dchem, rct, rcphot 
  use Zmet_ml  !> provides t, h2o etc from ESX variables

  implicit none
  private

  public :: ChemRun

  contains
  !-------------------------------------------------------------------------!
   subroutine ChemRun( dtchem, debug_level )
      real,    intent(in) :: dtchem
      integer, intent(in) :: debug_level
      integer :: i, k, nz, idj
      logical, save :: first_call = .true.
      real :: ppb, Jrate, ZenRad
      real :: fakenoon = 12*3600.0 ! for testing photorates only
      integer, save ::  Nout = 10  ! for debug output, Max NSPEC_TOT


      ppb = 1.0e9/Zmet(1)%M ! Output only


      if( first_call) then
        Nout = min(Nout, esx%nOutSpecs)

        call InitPhotol()
        call SetPhotolUsed()

        if( debug_level > 0 )  then
          write(*,"(a,2es12.2)") "CHDBG M, ppb ", Zmet(1)%M, ppb
          write(*,"(a,a7,20a8)") "CHDBG", "Time ", (species(i)%name, i=1, Nout)
          write(*,"(a,f6.0,20es8.1)") "CHDBG", 0.0, (xChem(i,1)*ppb, i=1,Nout)
        end if

      end if

      nz=esx%nz

      ZenRad = Loc%ZenDeg*DEG2RAD
      do i = 1, NPHOTOLRATES  !!! CRUDE. Needs also Z-profile
!PHUX        call setphotorates( photol_used(i), fakenoon, Jrate, debug_level)
        idj=photol_used(i)
        rcphot(idj,:) = GetPhotol(idj,ZenRad,debug_level )
        if( first_call .and. debug_level>0 ) &
          print "(a,f8.3,2i4,es14.4)", "PHOTOL ", Loc%ZenDeg, i, idj, rcphot(idj,1)
      end do

      call setchemrates( debug_level )

      if ( first_call ) then !.and. if ( debug_level > 0 ) then 
        write(*,*) "DEBUG RCT sizes: ", size(temp), size(rct(1,:)), size(rct(:,1))
        write(*,*) "DEBUG RCT k=1,10: "
        do k = 1, min( 10,  size(temp) )
          write (*,"(a,i3,2f6.1,3es10.3,2x,50es10.2)") "DEBUG RCT:", k, &
                TEMP(k), 1.0/TINV(k), M(k), O2(k), H2O(k), rct(1:4,k), rcphot(IDNO2,k)
        end do
      end if


    !> Do chemistry now, layer by layer for one dt_advec

      if(first_call) write(*, *) "ChemRun; nspec,nz=",  size(xChem,dim=1), size(xChem,dim=2)
      if( first_call) write(*,"(a,a6,20a8)") "CHDBG0", "Time", (trim(species(i)%name),i=1,Nout)
      if( first_call .or. debug_level > 1) then
         write(*,"(a,f6.0,20es8.1)") "CHDBG0", &
           (xChem( esx%OutSpecs(i)%int,1 )*ppb, i=1, Nout)
      end if

!print "(a,9es13.5)", "CHEMPRE", xChem(18,5), Dchem(18,5)
      do k = 1, Nz
          call chemsolve( k, dtchem, xChem(:,k),Dchem(:,k),debug_level )
      end do ! k
!print "(a,9es13.5)", "CHEMPOS", xChem(18,5), Dchem(18,5)

      if( debug_level > 1 ) &
        write(*,"(a,f6.0,20es8.1)") "CHDBG:", esx%Time, &
           (xChem( esx%OutSpecs(i)%int,1 )*ppb, i=1, Nout)

      first_call = .false.
     
   end subroutine ChemRun

end module esx_ChemRun
