module emep_Config_mod
  !---------------------------------------------------------------------
  ! A start of a more general module to store variables which
  ! can be easily changed though the namelist system. Will
  ! take over some of the job of the current ModelConstants;
  ! the latter stopped using just constants years ago ;-)
  !---------------------------------------------------------------------
  implicit none
  private

  type, private :: PBL_t
    real :: ZiMIN = 100.0                     ! minimum mixing height
    real :: ZiMAX = 3000.0                    ! maximum mixing height
    character(len=10) :: HmixMethod = "JcRb"  ! Method used for Hmix
      ! JcRb = Jericevic/Richardson number method
      ! "SbRb"= Seibert !"TIZi" = Original from Trond Iversen tiphysics
  end type PBL_t
  type(PBL_t), public, save :: PBL = PBL_t()
  

  type, private :: EmBio_t
    character(len=10) :: GlobBvocMethod = '-' ! can be MEGAN
    real :: IsopFac = 1.0                     ! for experiments
    real :: TerpFac = 1.0                     ! for experiments
  end type EmBio_t
  type(EmBio_t), public, save :: EmBio = EmBio_t()
  

end module emep_Config_mod
