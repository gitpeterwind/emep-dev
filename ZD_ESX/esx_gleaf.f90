!> MODULE esx_gleaf
!! Should call gleaf from relevant system. May need to re-structure
!! as part of DO3SE_ml?
module esx_gleaf

  use CheckStops, only: CheckStop_TF, CheckStop
  use esx_Variables, only: esx, Zveg, Loc, NZMAX => ESX_MAXNZ
  use esx_Zveg, only: Veg
  use esx_DO3SE, only: config_DO3SE, DO3SE_gleaf

  implicit none
  private

  public :: Set1Dgleaf

contains

  !> Sets 1-D leaf conductances, including gsto and gns components
  !! 

  subroutine Set1Dgleaf()
    integer :: nv, iz
    logical, save :: first_call = .true.

    nv = esx%nhVeg
    Zveg%gleaf = 0.0

    !print *, "GLEAFPARsun ", Loc%PARsun
    select case ( esx%gleaf_method )
     case ( 'gleaf0' ) 

        if ( Loc%PARsun > 1.0e-5 ) then
          Zveg(1:nv)%gleaf = Veg%gMax * ( 1.-exp(-Veg%alpha*Zveg(1:nv)%PARz ))
        end if

     case ( 'do3se' ) 

       if( first_call ) call config_DO3SE()
       call DO3SE_gleaf()                  !! Sets Zveg()%gsto
       Zveg(1:nv)%gleaf = Zveg(1:nv)%gsto  !! TMP

     case default
       call CheckStop("gleaf: unknown method: "//trim(esx%gsto_method))
     end select
     if ( esx%debug_gleaf > 0 ) then
       do iz = nv, 1, -1
           print *, "GLEAF Set1D iz PAR ("//trim(esx%gsto_method)//")", &
              iz, Zveg(iz)%PARz, Zveg(iz)%gleaf 
       end do
     end if
  end subroutine Set1Dgleaf

end module esx_gleaf
