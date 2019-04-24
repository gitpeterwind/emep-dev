!-------------------------------------------------------------------------!
!     Mapping of emissions from FINNv1 to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
      BiomassBurningMapping = "FINNv1.5toEmChem09"
  !FINNv1.5 uses 7 more GEOS-CHEM species and 2 more EMEP:
  integer, private, parameter :: &
    NBB_DEFS  = 15+7 & ! No mapping lines below
   ,NEMEPSPECS  = 13+2   ! No EMEP chemical mech specs used

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
  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  ",0.028 , 1.000, CO          ) &
   ,bbtype("NO  ",0.030 , 1.000, NO          ) &
   ,bbtype("NO2 ",0.046 , 1.000, NO2         ) &
   ,bbtype("SO2 ",0.064 , 1.000, SO2         ) &
   ,bbtype("NH3 ",0.017 , 1.000, NH3         ) &
   ,bbtype("ACET",0.058 , 1.000, C2H6        ) & !acetone 
   ,bbtype("ALD2",0.044 , 1.000, CH3CHO      ) &
   ,bbtype("ALK4",0.058 , 1.000, NC4H10      ) &
   ,bbtype("C2H6",0.030 , 1.000, C2H6        ) &
   ,bbtype("C3H8",0.044 , 0.700, NC4H10      ) & ! obs
   ,bbtype("CH2O",0.030 , 1.000, HCHO        ) &
   ,bbtype("MEK ",0.072 , 1.000, MEK         ) &
   ,bbtype("PRPE",0.042 , 1.000, C3H6        ) & 
   ,bbtype("PM25",1.0   , 1.000, PPM25_FIRE   ) & ! Tricky
  ! FINN v1.5, GEOS-CHEM 2015 changes: excludes 1 species:
  ! and needs one less EMEP species, C5H8
   ! ,bbtype("ISOP",0.068 , 1.000, C5H8        ) &
  ! FINN v1.5, GEOS-CHEM 2015 changes: added 8 species:
  ! and needs 3 more EMEP species = C2H4, OXYL, MGLYOX
   ,bbtype("C2H4",0.028 , 1.000, C2H4        ) & ! v1.5, new EMEP
   ,bbtype("C3H8",0.044 , 1.000, NC4H10      ) & ! v1.5
   ,bbtype("GLYC",0.060 , 1.000, CH3CHO      ) & ! v1.5 hydroxy-aceteldehyed
   ,bbtype("HAC", 0.060 , 1.000, C2H6        ) & ! v1.5 hydroxy-acetone. Query surrogate?
   ,bbtype("BENZ", 0.078 , 1.000, OXYL      ) & ! v1.5, new EMEP
   ,bbtype("TOLU", 0.092 , 1.000, OXYL      ) & ! v1.5
   ,bbtype("XYLE", 0.106 , 1.000, OXYL      ) & ! v1.5
   ,bbtype("MGLY", 0.072 , 1.000, MGLYOX      ) & ! v1.5, new EMEP
 /)
! Not used in EmChem09 standard
!   ,bbtype("OC  ",1.0   , 1.000, OC          ) &
!   ,bbtype("BC  ",1.0   , 1.000, BC          ) &
 !,bbtype("CH4  "      , 1.000, CH4         ) &
 ! ,bbtype("HCN "       , 1.000, HCN         ) &
