!-------------------------------------------------------------------------!
!     Mapping of emissions from FINNv1 to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
      BiomassBurningMapping = "FINNv1toEmChem09"
  integer, private, parameter :: &
    NBB_DEFS  = 15 & ! No mapping lines below
   ,NBBSPECS  = 15 & ! No  unique FINN used
   ,NCMSPECS  = 15   ! No EMEP chemical mech specs used

  !-----------------------------------------------------------------------!
  ! Column-1 gives GFED species ,
  ! Column-2 gives MW if needed to get to kg/day. Not needed for GFED
  !  (Some FINN emissions are given with mass basis, kg/day, so just
  !    set to 1.0. Others are given as mole/day, so multiply by MW/1000
  !    to get kg/day)
  ! Column-3 gives mass fraction of GFED emissions assign to EMEP species
  ! Column-4 gives EMEP species 
  !-----------------------------------------------------------------------!

  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  ",0.028 , 1.000, CO          ) &
   ,bbtype("NO  ",0.030 , 1.000, NO          ) &
   ,bbtype("NO2 ",0.046 , 1.000, NO2         ) &
   ,bbtype("SO2 ",0.064 , 1.000, SO2         ) &
   ,bbtype("NH3 ",0.017 , 1.000, NH3         ) &
   ,bbtype("ACET",0.058 , 1.000, C2H6        ) & ! CHECK acetone or acetylene?
   ,bbtype("ALD2",0.044 , 1.000, CH3CHO      ) &
   ,bbtype("ALK4",0.058 , 1.000, NC4H10      ) &
   ,bbtype("C2H6",0.030 , 1.000, C2H6        ) &
   ,bbtype("C3H8",0.044 , 0.700, NC4H10      ) & ! obs
   ,bbtype("CH2O",0.030 , 1.000, HCHO        ) &
   ,bbtype("ISOP",0.068 , 1.000, C5H8        ) &
   ,bbtype("MEK ",0.072 , 1.000, MEK         ) &
   ,bbtype("PRPE",0.042 , 1.000, C3H6        ) & !CHECK!
   ,bbtype("PM25",1.0   , 1.000, PPM25_FIRE   ) & ! Tricky
 /)
! Not used in EmChem09 standard
!   ,bbtype("OC  ",1.0   , 1.000, OC          ) &
!   ,bbtype("BC  ",1.0   , 1.000, BC          ) &
 !,bbtype("CH4  "      , 1.000, CH4         ) &
 ! ,bbtype("HCN "       , 1.000, HCN         ) &
