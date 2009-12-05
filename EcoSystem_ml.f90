module EcoSystem_ml

 use CheckStop_ml,       only : CheckStop, StopAll
 use GridValues_ml,      only : debug_proc, debug_li, debug_lj
 use LandDefs_ml,        only : LandDefs, LandType
!CIRCULAR SOMEHERE use Landuse_ml !!!,         only : LandCover
 use ModelConstants_ml , only : MasterProc, DEBUG_ECOSYSTEMS, NLANDUSEMAX
 use OwnDataTypes_ml,    only : Deriv, print_deriv_type &
                                ,TXTLEN_DERIV, TXTLEN_SHORT
 use Par_ml,             only : li0, lj0, li1, lj1, MAXLIMAX, MAXLJMAX
 implicit none
 private

    public :: Init_EcoSystems

  ! depositions are calculated to the following landuse classes, where
  ! e.g. conif may include both temperate and Medit. forests

  ! Water_D
  ! *** Note *** Water_D is introduced for some NEU work, with direct 
  ! deposition to the water surface. This is not to be used for IAM, 
  ! since CCE want to have deposition to the watershed, which means 
  ! the grid in practice.

   integer, parameter, public :: NDEF_ECOSYSTEMS = 6
   character(len=8), public, dimension(NDEF_ECOSYSTEMS), parameter :: &
    DEF_ECOSYSTEMS = (/ "Grid    " &
                      , "Conif   " &
                      , "Decid   " &
                      , "Crops   " &
                      , "Seminat " &
                      , "Water_D " /)
   type(Deriv), public, dimension(NDEF_ECOSYSTEMS), save :: DepEcoSystem

   logical,  public, dimension(NDEF_ECOSYSTEMS,NLANDUSEMAX), &
            save :: Is_EcoSystem

   real, public, dimension(NDEF_ECOSYSTEMS,MAXLIMAX,MAXLJMAX), &
            save :: EcoSystemFrac

   integer, public, parameter :: FULL_GRID=1
   integer, private, parameter :: &
     CONIF=2, DECID=3, CROP=4, SEMINAT=5, WATER_D=6 ! try to skip

contains
 !<---------------------------------------------------------------------------
  subroutine Init_EcoSystems()

    real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
    character(len=TXTLEN_DERIV) :: name
    character(len=TXTLEN_SHORT)  :: unit
    integer :: i,j,ilc,lc,nlc, iEco
    logical, parameter :: T = .true., F = .false. ! shorthands only
    logical :: debug_flag
    real :: coverage

    do iEco = 1, NDEF_ECOSYSTEMS

       if( MasterProc ) &
           write(*,*) "======== DEP REECEIVERS " // DEF_ECOSYSTEMS(iEco)

       name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_Frac"
       unit = "Fraction"
       if(iEco==FULL_GRID) then
          name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_km2"
          unit = "km2"
       end if

       ! Set defaults
        ! dep_type( name, LC, index, f2d, class, label, txt, scale, atw, units )
        !            x     d      d    d   a10    a10   a10     f    i    a10

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,  LC, XYLC, scale, dt_scale, avg? rho Inst Yr Mn Day atw
        DepEcoSystem(iEco) = Deriv(  &
               name, "EcoFrac", "Area",DEF_ECOSYSTEMS(iEco) , unit, &
                  iEco, -99, iEco,-99.9,  1.0, 0.0,  F,   F , F ,T ,F ,F, -99.9 )


        if(DEBUG_ECOSYSTEMS .and. MasterProc) &
             call print_deriv_type( DepEcoSystem(iEco) )

    end do

 !  Define which landcovers belong to which ecosystem

        Is_EcoSystem(FULL_GRID,:)    =  .true.
        Is_EcoSystem(CONIF,:)   =  LandType(:)%is_conif
        Is_EcoSystem(DECID,:)   =  LandType(:)%is_decid
        Is_EcoSystem(CROP,:)    =  LandType(:)%is_crop
        Is_EcoSystem(SEMINAT,:) =  LandType(:)%is_seminat
        Is_EcoSystem(WATER_D,:) =  LandType(:)%is_water

    EcoSystemFrac(:,:,:) = 0.0

  end subroutine Init_EcoSystems

end module EcoSystem_ml
