!-------------------------------------------------------------------------!
!     Mapping of emissions from FINNv1 to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
      BiomassBurningMapping = "DUMMY"
  integer, private, parameter :: &
    NBB_DEFS  = 1 & ! No mapping lines below
   ,NBBSPECS  = 1 & ! No  unique FINN used
   ,NCMSPECS  = 1   ! No EMEP chemical mech specs used

  !-----------------------------------------------------------------------!
  ! Column-1 gives FINNv1 species ,
  ! Column-2 gives MW if needed to get to kg/day. Not needed for FINNv1
  !  (Some FINN emissions are given with mass basis, kg/day, so just
  !    set to 1.0. Others are given as mole/day, so multiply by MW/1000
  !    to get kg/day)
  ! Column-3 gives mass fraction of FINNv1 emissions assign to EMEP species
  ! Column-4 gives EMEP species 
  !-----------------------------------------------------------------------!

  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  ",0.028 , 1.000, CO          )  &
 /)
! Not used in EmChem09 standard
!   ,bbtype("OC  ",1.0   , 1.000, OC          ) &
!   ,bbtype("BC  ",1.0   , 1.000, BC          ) &
 !,bbtype("CH4  "      , 1.000, CH4         ) &
 ! ,bbtype("HCN "       , 1.000, HCN         ) &
