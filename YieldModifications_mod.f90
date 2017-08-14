! <YieldModifications_mod.f90 - part of the EMEP MSC-W Chemical transport Model>
!_____________________________________________________________________________!
! Added March 2017
! Can be empty!
! For placement of miscellaneous coefficients and routines which can be
! called from Solver_ml. Currently allows access to JPC-based yields
! D. Simpson March 2017
!_____________________________________________________________________________!
!
 module YieldModifications_mod

  use CheckStop_ml,       only : StopAll
  use ChemFields_ml             ! => NSPEC_TOT, O3, NO2, etc.
  use ChemSpecs                  ! => NSPEC_TOT, O3, NO2, etc.
  use emep_Config_mod,    only : YieldModifications
  use ModelConstants_ml,  only : MasterProc, DebugCell, DEBUG, USES! JPCIsoYield
  use NumberConstants,    only : UNDEF_R, UNDEF_I
  use SmallUtils_ml,      only : find_index, trims

  implicit none
  private

  public  :: doYieldModifications

  ! Coded up some JPC experiments for now
  ! Allows JPC-FY-acid or JPC-VY-neutral or similar

  private :: init_JPCyields
  private :: update_JPCyields

  logical, public, save :: YieldModificationsInUse = .false.

   real, public, save :: & ! CRUDE, but fix later
      YA0APINOH   &  ! Limit of 17% mass at low Isop
     ,YG_APINOH, YG_APINO3, YG_APINNO3, YA_APINOH, YA_APINO3, YA_APINNO3   &
     ,YG_BPINOH, YG_BPINO3, YG_BPINNO3, YA_BPINOH, YA_BPINO3, YA_BPINNO3   &
     ,YG_MTOH,   YG_MTO3,   YG_MTNO3,   YA_MTOH,   YA_MTO3,   YA_MTNO3  &
     ,YA_ISOPOH, YA_ISOPNO3 

   real, private, save    :: IsoOHyield = UNDEF_R  ! 1 or 4%  for JPC


 contains

  ! --------------------------------------------------------------------------
  !> SUBROUTINE doYieldModifications
  !! Called once on Solver first_call (which sets YieldModificationsInUse, and some
  !! constants yield values), at start of each set of chemical iterations
  !! for each grid-cell (to reset yields for this cell), then after each iteration 
  !! to update the yields based.

  subroutine doYieldModifications(txt)
     character(len=*), intent(in) :: txt

     if ( YieldModifications(1:3) == 'JPC' ) then

       call init_JPCyields(txt)

       if ( YieldModifications(1:6) == 'JPC-VY' .and. txt == 'run' )  then
         call update_JPCyields()  ! OH impact
       end if

     end if
  end subroutine doYieldModifications

  ! --------------------------------------------------------------------------
  !> SUBROUTINE init_JPCyields
  !!
  subroutine init_JPCyields(txt)
  !>--------------------------------------------------------------------------
   character(len=*), intent(in) :: txt
   character(len=*), parameter :: dtxt='initJPC:'
   logical, save :: first_call = .true.

 ! WARNING hard-coded MWs here to allow compilation with EmChem. Will fix later
 ! For isoprene SOA product OM/OC ratio set to 2.0 (as in Bergstrom et al., 
 ! 2012 - based on Chhabra et al., 2010), which gives SOA MW=136, ISOP MW = 68
 ! For MT SOA we use MW=240

   real, parameter :: MWratioISOP = 68.0/136.0  ! = 0.5
   real, parameter :: MWratioMT = 136.0/240.0   ! = 0.5667


   if ( first_call ) then

     YieldModificationsInUse = .true.
     if( MasterProc ) write(*,*) dtxt//'YieldModifications:'//&
           trim(YieldModifications)

     if ( index( YieldModifications, 'acid') > 0 ) then
         IsoOHyield = 0.04                   ! default is acid, high yield ?
     else ! neutral
         IsoOHyield = 0.01 
     end if

   ! Isoprene

      YA_ISOPOH  = IsoOHyield * MWratioISOP

   ! Assume 4% (mass-based) BSOA yield for isoprene + NO3 reaction - based
   ! on lowest yield observed by Ng et al. 2008

      YA_ISOPNO3 = 0.04 * MWratioISOP

   ! Monoterpenes. Yields depend on dC5H8/dMT if JPC-VY

     ! OH: 
      YA0APINOH = 0.17 * MWratioMT !species(APINENE)%molwt/species(APINOH_BSOA_NV)%molwt

     ! O3, fixed yields
      YA_APINO3 = 0.15 * MWratioMT !species(APINENE)%molwt/species(APINO3_BSOA_NV)%molwt
      YG_APINO3 = 1 - YA_APINO3

      YA_BPINO3 = YA_APINO3
      YG_BPINO3 = YG_APINO3
      YA_MTO3   = YA_APINO3
      YG_MTO3   = 1- YA_APINO3

     ! NO3, fixed yields
      YA_BPINNO3 = 0.3  * MWratioMT
      YG_BPINNO3 = 1 - YA_BPINNO3

      YA_MTNO3   =  0.3 * MWratioMT
      YG_MTNO3   =  1- YA_APINNO3

      first_call = .false.

      if ( MasterProc ) write(*,*) dtxt//trim(txt), IsoOHYield, YA0APINOH
    end if

    ! Needed at start of each chem timestep (may have changed if JPC-VY used)
    YA_APINOH = YA0APINOH
    YG_APINOH = 1 - YA_APINOH

    YA_BPINOH = YA_APINOH
    YG_BPINOH = YG_APINOH
    YA_MTOH = YA_APINOH
    YG_MTOH = 1- YA_APINOH

    if ( DEBUG%RUNCHEM .and. DebugCell ) then
       write(*,"(3a,3es12.3)") dtxt//trim(txt), &
           DEBUG%datetxt, "====", YA_MTOH, YA_MTO3,YA_MTNO3
    end if
    
  end subroutine init_JPCyields

  ! --------------------------------------------------------------------------
  !> SUBROUTINE update_JPCyields
  !  Only needed for JPC VY runs, and only affects OH yields

   subroutine update_JPCyields()
     character(len=*), parameter :: dtxt='updateJPCY:'
     real, parameter :: Y1=0.48, YD = 1-Y1
     integer, save :: i_iso=UNDEF_I, i_mt=UNDEF_I ! indices of OHLOSS
     logical, save :: first_call = .true.
     real          ::  diso, dmt, ratio, yield

     if ( first_call ) then
        i_iso = find_index('OHLOSS_ISO',species(:)%name)
        i_mt  = find_index('OHLOSS_MT',species(:)%name)
        if(MasterProc) write(*,*) dtxt//' i   ', i_iso, i_mt
        if ( i_iso < 1 .or. i_mt < 1 ) then
           print *, dtxt//"MISSING OHLOSS for JPC: ", i_iso, i_mt
           call StopAll(dtxt//'JPC OHLOSS ERR')
       end if
       first_call = .false.
     end if

     ! nb 1.0e5 OH * 10 ppt C5H8 * rc  -> P(OHLOSS)~2500 

     diso=xnew(i_iso)

     if ( diso < 1.0 ) return
      
     dmt =xnew(i_mt)
     ratio = diso /(1.0+diso + dmt )

     if ( DEBUG%RUNCHEM .and. DebugCell ) then
        write(*,"(4a,4es12.3)") dtxt, DEBUG%datetxt, "====", &
             dtxt, diso,dmt, ratio, xnew(C5H8)
     end if

     if ( ratio < 0.01 )  then    ! exact to 3 sig figs
         yield  = 1 - ratio
     else if ( ratio < 5 )  then    ! exp(-1.53*5) = 2.5e-4 anyway: y=0.4802
         !y =          0.48 + 0.52 * exp(-1.53 * ratio )
         yield  =        Y1   + YD * exp(-1.53 * ratio )
     else 
         yield = Y1
     end if

    YA_APINOH = YA0APINOH * yield
    YG_APINOH = 1 - YA_APINOH

    ! We assume same OH behaviour for other MT:
    YA_BPINOH = YA_APINOH
    YG_BPINOH = 1-YG_BPINOH

    YA_MTOH = YA_APINOH
    YG_MTOH = 1- YA_MTOH

    if ( DEBUG%RUNCHEM .and. DebugCell ) then
      write(*,"(2a,9es12.3)") dtxt//' diso dmt rat y ', DEBUG%datetxt, diso, dmt, ratio, yield
      if ( ratio > 0.01 .and. ratio < 0.5  ) then
            write(*,"(3a,9es12.3)") dtxt//"JPC5YIELD ",dtxt, DEBUG%datetxt, &
              diso, dmt, ratio, yield, YA_APINOH
      end if
    end if

  end subroutine update_JPCyields

 end module YieldModifications_mod
