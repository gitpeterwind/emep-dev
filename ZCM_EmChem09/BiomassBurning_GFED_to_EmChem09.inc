
!-------------------------------------------------------------------------!
!     Mapping of emissions from GFED to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!
  character(len=*), public, parameter :: &
     BiomassBurningMapping = "GFEDtoEmChem09"
  integer, private, parameter :: &
    NBB_DEFS  = 6 & ! No mapping lines below
   ,NBBSPECS  = 5 & ! No  unique GFED used
   ,NCMSPECS  = 6 & ! No EMEP chemical mech specs used
   ,MAXBB_SPLIT= 2   ! Max No EMEP components per GFED, e.g for VOC
 
  !-----------------------------------------------------------------------!
  ! Column-1 gives GFED species , 
  ! Column-2 gives MW if needed to get to kg/day. Not needed for GFED
  !  (GFED emissions are given with mass basis, so just set to 1.0)
  ! Column-3 gives mass fraction of GFED emissions assign to EMEP species
  ! Column-4 gives EMEP species
  !-----------------------------------------------------------------------!

  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  " , 1.0 , 1.000, CO          ) &
   ,bbtype("PM25" , 1.0 , 1.000, PPM25_FIRE   ) & ! Tricky
   ,bbtype("NMHC" , 1.0 , 0.700, NC4H10       ) & ! TMP!!!
   ,bbtype("NMHC" , 1.0 , 0.300, C2H6         ) & ! TMP!!!
   ,bbtype("NOx " , 1.0 , 1.000, NO          ) &
   ,bbtype("SO2 " , 1.0 , 1.000, SO2         ) &
 /)
! Not used in EmChem09 standard
 !,bbtype("TPM         , 1.000, CH4         ) &
!   ,bbtype("OC  "       , 1.000, OC          ) &
!   ,typ_sfi("BC  "       , 1.000, BC          ) &
