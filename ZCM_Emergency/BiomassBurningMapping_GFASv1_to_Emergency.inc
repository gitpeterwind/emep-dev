!-------------------------------------------------------------------------!
!     Mapping of emissions from GFASv1 to EMEP Emergency species
!     Included by  ForestFire_ml
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
    BiomassBurningMapping = "GFASv1toEmergency"
  integer, private, parameter :: &
    NBB_DEFS   = 0, & ! No. mapping lines below
    NEMEPSPECS = 0    ! No. EMEP chemical mech specs used

!-----------------------------------------------------------------------!
! Column-1 gives GFASv1 species,
! Column-2 gives MW if needed to get to kg/day. Not needed for GFASv1
!  (GFAS units are [kg/m2/s]. No unit conversion is needed.)
! Column-3 gives mass fraction of GFASv1 emissions assign to EMEP species
! Column-4 gives EMEP species
!-----------------------------------------------------------------------!
  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs
