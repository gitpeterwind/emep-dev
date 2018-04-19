
!-------------------------------------------------------------------------!
!     Mapping of emissions from GFED to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!
  character(len=*), public, parameter :: &
     BiomassBurningMapping = "GFEDtoEmChem09ecom"
  integer, private, parameter :: &
    NBB_DEFS  = 15 & ! No mapping lines below
   ,NEMEPSPECS  = 13 !& ! No EMEP chemical mech specs used
!   ,MAXBB_SPLIT= 7   ! Max No EMEP components per GFED, e.g for VOC
 
  !-----------------------------------------------------------------------!
  ! Column-1 gives GFED species , 
  ! Column-2 gives MW if needed to get to kg/day. Not needed for GFED
  !  (GFED emissions are given with mass basis, so just set to 1.0 except 
  !   for OC, which is given in units of C and needs to be scaled from OC 
  !   to OM)
  ! Column-3 gives mass fraction of GFED emissions assign to EMEP species
  ! Column-4 gives EMEP species
  !-----------------------------------------------------------------------!

  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  " , 1.0 , 1.00, CO          ) &
   ,bbtype("PM25" , 1.0 , 1.00, FFIRE_REMPPM25 ) & ! Tricky, -(OM+BC)
   ,bbtype("NMHC" , 1.0 , 0.20, C2H6         ) & 
   ,bbtype("NMHC" , 1.0 , 0.10, NC4H10       ) & 
   ,bbtype("NMHC" , 1.0 , 0.20, C2H4         ) & 
   ,bbtype("NMHC" , 1.0 , 0.20, C3H6         ) & 
   ,bbtype("NMHC" , 1.0 , 0.20, OXYL         ) & 
   ,bbtype("NMHC" , 1.0 , 0.05, HCHO         ) & 
   ,bbtype("NMHC" , 1.0 , 0.05, CH3CHO       ) & 
   ,bbtype("NOx " , 1.0 , 1.000, NO          ) &
   ,bbtype("SO2 " , 1.0 , 1.000, SO2         ) &
   ,bbtype("OC  " , 1.7 , 1.000, FFIRE_OM    ) &
   ,bbtype("BC  " , 1.0 , 1.000, FFIRE_BC    ) &
   ,bbtype("OC  " ,-1.7 , 1.000, FFIRE_REMPPM25    ) &  ! Subtract from PM25 used above
   ,bbtype("BC  " ,-1.0 , 1.000, FFIRE_REMPPM25    ) &  ! Subtract from PM25 used above
 /)
! Not used in EmChem09 standard
 !,bbtype("TPM         , 1.000, CH4         ) &
