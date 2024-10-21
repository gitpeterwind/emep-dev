!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************!

module PBAP_mod
  !/-- Module to deal with primary biological aerosol particles (PBAPs).
  !    DUMMY! Will be updated and released in due course
  !
  !    Gunnar Felix Lange 2024
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use Config_module, only : USES

  implicit none

  public ::  init_PBAPs,set_PBAPs

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    if (USES%FUNGAL_SPORES) then
      call StopAll('ERROR: You have set USE%FUNGAL_SPORES = T, but this is not currently implemented.')
    end if
    if (USES%BACTERIA) then
      call StopAll('ERROR: You have set USE%BACTERIA = T, but this is not currently implemented.')
    end if
    if (USES%MARINE_OA) then
      call StopAll('ERROR: You have set USE%MARINE_OA = T, but this is not currently implemented.')
    end if
  end subroutine init_PBAPs


  subroutine set_PBAPs(i,j)
    integer, intent(in) ::  i,j
    call StopAll('ERROR: You are trying to call Primary Biological Aerosol Particles (PBAPs),but they are not currently implemented.')

  end subroutine set_PBAPs

end module PBAP_mod
