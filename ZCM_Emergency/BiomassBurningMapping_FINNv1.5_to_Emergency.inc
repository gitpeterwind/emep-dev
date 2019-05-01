!-------------------------------------------------------------------------!
!     Mapping of emissions from FINNv1 to EMEP Emergency species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
      BiomassBurningMapping = "FINNv1.5toEmergency"
  integer, private, parameter :: &
    NBB_DEFS  = 0 & ! No. mapping lines below
   ,NEMEPSPECS  = 0   ! No. EMEP chemical mech specs used

  !-----------------------------------------------------------------------!
  ! Column-1 gives FINNv1 species ,
  ! Column-2 gives MW if needed to get to kg/day. Not needed for FINNv1
  !  (Some FINN emissions are given with mass basis, kg/day, so just
  !    set to 1.0. Others are given as mole/day, so multiply by MW/1000
  !    to get kg/day)
  ! Column-3 gives mass fraction of FINNv1 emissions assign to EMEP species
  ! Column-4 gives EMEP species 
  !-----------------------------------------------------------------------!

  !     BBname    unitsfac frac  emep:
  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs
