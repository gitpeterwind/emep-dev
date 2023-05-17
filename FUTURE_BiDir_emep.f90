! <BiDir_emep.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33BiDir2>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
! CODE is merge of original BiDir_emep, DryDep_module from Roy/Dave and later
! work from Sebastiaan. 
!*****************************************************************************!

module BiDir_emep
  use AllocInits,         only: AllocInit
  use BiDir_module,       only: BiDir, BiDir_2d, &
             BiDirXconcs, BiDirResistances, BiDirFluxes, &
             BiDir_errmsg, BiDirXwaterEuro, BiDirXwaterOrig
  use Biogenics_mod,        only: SoilNH3  ! for BiDir
  use CheckStop_mod
  use ChemSpecs_mod,      only: NH3
  use Config_module,      only: KMAX_MID,MasterProc, USES, IOU_INST, &
                                step_main, iyr_trend ! QUERY yrtrend or meteo yr?
  use Debug_module,       only: DEBUG
  use DerivedFields_mod,  only: d_2d
  use GridValues_mod,     only: debug_proc, debug_li,debug_lj, glon, glat &
                                       ,i_fdom, j_fdom ! DS added
  use Io_Progs_mod,       only: datewrite
  use LandDefs_mod,       only: LandDefs  !DS added
  use LocalVariables_mod, only: Grid, L     ! Meteorology as e.g. Grid%t2C
  use NetCDF_mod,         only: ReadField_CDF, printCDF
  use PhysicalConstants_mod, only: AVOG
  use Rsurface_mod,       only: Rinc !DS added
  use ZchemData_mod,      only: xn_2d
  use PAR_mod,            only: LIMAX, LJMAX, me
  use SmallUtils_mod,     only: key2str
  use SubMet_mod,         only: Sub  ! DS added
  use TimeDate_mod,       only: date, current_date, print_date

  implicit none
  private

  private :: BiDir_InitFields  ! Allocates and sets long-term NH3,
  public :: BiDir_ijInit     ! updates monthly fields, then calls BiDirXconcs
  public :: BiDir_ijRGs      ! calls to Resistances
  public :: BiDir_ijFluxes   ! calls to BiDirFluxes
  public :: BiDir_ijFinish   ! reset gradient_fac, saves concs, and update SoilNH3
  public :: BiDir_Derived    ! sets d_2d for e.g. ugXTOT
  ! If USES%BiDirEuroXwater:
   public :: load_input_data_BiDir ! Allocates and sets long-term NH3, updates monthly fields
   public :: set_BiDirSea
   private :: fixdates
  integer, public :: old_month
  
! moved here from DryDep
  real, private, save :: &
    wet, dry        & ! Fractions
   ,Xtot, Xemis, GextBD, GigsBD, ug_old   &
   ,ugNH3_zmid, ugNH3_3m_LC, ugXH3_3m_LC, ugNH3_3m_grid, ugXH3_3m_grid &
   ,sstK, ssNH4, sspH, ssS, ss, Xwater
  real, parameter ::  xn2ugNH3 =  17/AVOG*1.0e6*1.0e6  ! molec/cm3  -> ug/m3
  real, save :: oldRsur
  real, allocatable, dimension(:,:) :: tmpNHxDep !BIDIR for OrigXwater

  logical, private, save :: dbgBD

contains
!----------------------------------------------------------------------------
!========== INIT ON FIRST CALL =======================================
  subroutine BiDir_InitFields()
     integer :: i,j ! for debug
     character(len=*), parameter :: dtxt='BiDir:Fields:'
     
     if ( MasterProc) write(*,*) dtxt//'called DATE',step_main, iyr_trend, print_date()
     
     allocate (BiDir_2d%NH3aLT(LIMAX,LJMAX))
     allocate (BiDir_2d%SO2aLT(LIMAX,LJMAX))!Hazelhos 15-11-2019: added this line

     allocate (BiDir_2d%aSN(LIMAX,LJMAX))
     allocate (BiDir_2d%Xtot(LIMAX,LJMAX))
     allocate (BiDir_2d%NH3inst(LIMAX,LJMAX))
     allocate (BiDir_2d%NHxEmis(LIMAX,LJMAX))
     allocate (BiDir_2d%NH3_3m(LIMAX,LJMAX))
     allocate (BiDir_2d%XH3_3m(LIMAX,LJMAX))
     
     if (BiDir%OrigXwater) then
       allocate (BiDir_2d%NHxDep(LIMAX,LJMAX)) !DSMAY23
       allocate (tmpNHxDep(LIMAX,LJMAX)) !DSMAY23
     else if (BiDir%EuroXwater) then
       allocate (BiDir_2d%sea_nh4(LIMAX,LJMAX))
       allocate (BiDir_2d%sea_ph(LIMAX,LJMAX))
       allocate (BiDir_2d%sea_gridS(LIMAX,LJMAX))
       allocate (BiDir_2d%sea_gridT(LIMAX,LJMAX))
     end if
     
     BiDir_2d%NH3inst = -999.  ! Initial value as flag  0.0
     old_month = -999
     
    !Hazelhos 15-11-2019:
    !A standard simulation, specified by BiDir%InputFile in the namelist file, is used for reading SO2 and NH3 
    !values. This is a temporary solution, as we want to obtain SURF_ug_SO2 and SURF_ug_NH3 from the global
    !model for the first domain, and use for the other domains the values from domain 1
    
    
     ! 1) Annual NH3, ug/m3
     call ReadField_CDF(BiDir%InputFile, 'SURF_ug_NH3',&
        BiDir_2d%NH3aLT, nstart=1,interpol='zero_order',needed=.true.,&
        UnDef=0.0, debug_flag=.false.)

     call printCDF('BIDIR_NH3',BiDir_2d%NH3aLT,'ugm3')

     !     2) Annual SO2, ug/m3 !Hazelhos 15-11-2019: added
     call ReadField_CDF(BiDir%InputFile, 'SURF_ug_SO2',&
        BiDir_2d%SO2aLT, nstart=1,interpol='zero_order',needed=.true.,&
        UnDef=0.0, debug_flag=.false.)

     call printCDF('BIDIR_SO2',BiDir_2d%SO2aLT,'ugm3')
     
     
     !Hazelhos 15-11-2019: Outdated. is not used. Can be removed?
     !DS: recoded with USES%BiDirOrigWater, as an option for testing.
     !Longer term solution still unclear.
     
     !  2) Annual DDEP NHx, mgN/m2
     ! WHY? call ReadField_CDF(BiDir%InputFile,'DDEP_RDN_m2Water_D', BiDir_2d%NHxDep,&
     if (BiDir%OrigXwater) then
       call ReadField_CDF(BiDir%InputFile,'WDEP_RDN', tmpNHxDep,& 
           nstart=1,interpol='zero_order',needed=.true.,UnDef=0.0, &
             debug_flag=.false.)
       call ReadField_CDF(BiDir%InputFile,'DDEP_RDN_m2Grid', BiDir_2d%NHxDep,& 
           nstart=1,interpol='zero_order',needed=.true.,UnDef=0.0, &
             debug_flag=.false.)
       BiDir_2d%NHxDep = BiDir_2d%NHxDep + tmpNHxDep

       ! from mgN/m2 to kg(NH4???)/ha:

       BiDir_2d%NHxDep = BiDir_2d%NHxDep * 1.0e-6*1.0e4*18.0/14  
     end if ! USES%BiDirOrigWater
     
     !Hazelhos 20191115: calculate aSN based on SO2 and NH3 (convert
     !  ug m-3 to mol m-3), rather than reading it directly from a global file
     BiDir_2d%aSN = (BiDir_2d%SO2aLT / 64.066) / (BiDir_2d%NH3aLT / 17.031)
     call printCDF('BIDIR_aSN',BiDir_2d%aSN,'ratio')
     
     if( debug_proc ) then
       i=debug_li
       j=debug_lj
       write(*, "(a,2i4,2f8.2)") dtxt//' coords/lonlat', &
          i_fdom(i),j_fdom(j), glon(i,j), glat(i,j)
       write(*, "(a,3a12)") dtxt//"    ",'ijVal','min','max'
       write(*, "(a,3f12.3)") dtxt//"NH3LT ",BiDir_2d%NH3aLT(i, j),&
            minval(BiDir_2d%NH3aLT(:,:)), maxval(BiDir_2d%NH3aLT(:,:))
       write(*, "(a,3f12.3)") dtxt//"SO2LT ",BiDir_2d%SO2aLT(i, j),&
            minval(BiDir_2d%SO2aLT(:,:)), maxval(BiDir_2d%SO2aLT(:,:))
       write(*, "(a,3f12.3)") dtxt//"aSN ",BiDir_2d%aSN(i,j),&
            minval(BiDir_2d%aSN(:,:)), maxval(BiDir_2d%aSN(:,:))
      end if

  end subroutine BiDir_InitFields
!----------------------------------------------------------------------------
  subroutine BiDir_ijInit(i,j,NH3_ix)
     integer, intent(in) :: i,j,NH3_ix
     character(len=20) :: BiDirMsg  ! issues warnings from BiDir if needed 
     character(len=*), parameter :: dtxt='BiDir:Init:'
     logical, save :: first_call = .true.
     real :: sstK, ssNH4, sspH, S

       if ( first_call ) then

         call BiDir_InitFields()

         first_call = .false.

       end if

       dbgBD = ( DEBUG%BIDIR .and. debug_proc .and.  i == debug_li .and.&
         j == debug_lj .and.  mod(current_date%hour,6) ==  0 )   ! every 6 hours

       ugNH3_zmid = xn_2d(NH3_ix,KMAX_MID) * xn2ugNH3  ! conc. at grid centre, ug/m3

       if( BiDir_2d%NH3inst(i,j) < 0 ) BiDir_2d%NH3inst(i,j) = ugNH3_zmid ! Initial

       ! Get Xwater, Xgs, Xext:
       Xwater = -999.0

       !if( USES%BiDirEuroXwater .and. L%is_water ) then !FIXME is_water???
       if( BiDir%EuroXwater ) then
          !Hazelhos 20191114 inserted load_input_data_BiDir to check if current
          ! month has changed and Sea data needs to be updated (load sea data)

          call load_input_data_BiDir()
          call set_BiDirSea(i,j,ssNH4,sspH,sstK,ssS)
         
          if ( dbgBD ) write(*,'(a, 2es12.3, L2)') "sstK ", &
                  sstK, Xwater, dbgBD

          !if( sstK < 500.0  ) & !Hazelhos 20-03-2020: removed if-statement.
          ! We want Xwaters to be calculated also for non-CMEM
          !cells. Now, Xwater is calculated for every cell,
          !including non-water cells. It doesn't have an effect
          !on non-water cells, since it is overwritten in BiDirXconcs.

          Xwater=BiDirXwaterEuro(i,j,Grid%sst,ssNH4,sspH,sstK,ssS,&
            Sub(14)%is_water,dbgBD,'nwpSST') !or monthlySST, 
             ! Hazelhos 21-02-2020: added Sub(14)%is_water boolean and dbgBD.
             ! 14 = water LU. Should make this more general?

          if ( dbgBD ) write(*,'(a,9es12.3)') "BDXwater ", ssNH4,sspH,&
            sstK,ssS, Grid%sst, Xwater

       else if ( BiDir%OrigXwater ) then !DSMAY23

         Xwater=BiDirXwaterOrig(i,j,Grid%sst,BiDir_2d%NHxDep(i,j))

       end if !USES%BiDirEuroXwater

       if ( dbgBD ) then
         write(*,'(a, 3L2,e12.3)') dtxt//print_date()//&
         " L%is_water, G%is_mainlysea, Sub(14)%is_water,Xwater "& 
          ,L%is_water, Grid%is_mainlysea, Sub(14)%is_water, Xwater
       end if

      ! ========================================= BIDIR =============!

       call BiDirXconcs(Grid%t2C, Grid%sst, Xwater, BiDir_2d%aSN(i,j),  &
           BiDir_2d%NH3inst(i,j), &  ! 3m conc (ug/m3) from last time-step
           BiDir_2d%NH3aLT(i,j), dbg=dbgBD)

       call CheckStop(BiDir_errmsg/='ok',BiDir_errmsg)
       ! ========================================= BIDIR =============!

       !Hazelhos 15-11-2019: removed BiDir_NHxDep(i,j), 
       !Hazelhos 07-02-2020:
       !  added water boolean Grid%is_mainlysea (is valid when fraction
       ! water > 0.5). May not be a sustainable solution, since we want to
       ! introduce river emissions and IJsselmeer as well. Not sure if
       ! IJsselmeer qualifies as Grid%is_mainlysea
       ! consider moving this function call for BiDirXconcs to below for a
       ! better solution, when we have a value for L%is_water.

       !Hazelhos 21-02-2020: changed Grid%is_mainlysea with Sub(14)%is_water.
       ! 14 is the land use category for water.
       !Hazelhos 20-03-2020: 
       !   Removed Sub(14)%is_water as it is no longer necessary.
       !DS if ( dbgBD ) write(*,'(a, 1es12.3)') "Xwater After Xconcs: ",  Xwater

      ! Initialise before LC loop:
       BiDir_2d%NHxEmis(i,j) = 0.0
       BiDir_2d%Xtot(i,j)    = 0.0
       BiDir_2d%NH3_3m(i,j)    = 0.0
       BiDir_2d%XH3_3m(i,j)    = 0.0
       ugNH3_3m_grid = 0.0 ! from zmid down
       ugXH3_3m_grid = 0.0 ! from Xtot upwards

  end subroutine BiDir_ijInit
!----------------------------------------------------------------------------
  subroutine BiDir_ijRGs(ncall,iL,Rsur,Gsto)
     integer, intent(in):: ncall, iL
     real, intent(out) :: Rsur
     real, intent(in)  :: Gsto
     character(len=*), parameter :: dtxt='BiDir_getRes:'
     character(len=30) :: otxt

      if ( dbgBD ) then ! Hazelhos prints (DS converted to write)
        oldRsur=Rsur
        write(otxt,"(a,2i3)") dtxt// 'inputs ', ncall, iL
        write(*,'(a,2L2,99(a,g10.3))') otxt// 'ice water:',  &
          Grid%is_frozen,L%is_water, ' L%rh', L%rh, ' G%rh2m',  Grid%rh2m
        write(*,'(a,99(a,g10.3))') otxt,  &
          ' Rinc', Rinc, ' Gsto', Gsto, 'oldRs',oldRsur
      end if

       !Hazelhos, autumn 2019: changed Grid%rh2m to L%rh, because Grid%rh2m
       ! gives fractions of order E-5 whereas L%rh gives more realistic
       ! fractions ~0.7

        call BiDirResistances(L%SAI,L%rh,Grid%is_frozen,L%is_water,&
          Rinc,Gsto, &                      !:inputs
            GextBD,GigsBD,Rsur,dbgBD)  !:outputs

        if ( dbgBD ) then
          write(otxt,"(a,i3)") dtxt// 'outputs ', ncall
          write(*,'(a, 3es9.2)') otxt// &
         ' GextBD, GigsBD, Rsur', GextBD, GigsBD, Rsur
        end if

  end subroutine BiDir_ijRGs

!----------------------------------------------------------------------------
  subroutine BiDir_ijFluxes(i,j,iL,Vg_ref,Vg_eff,Rb,Rsur,Gsto)
     integer, intent(in) :: i,j,iL
     real, intent(in) :: Vg_ref,Vg_eff,Rb,Rsur,Gsto
     real :: Ra_diff
     character(len=*), parameter :: dtxt='BiDirFlux'

     if ( dbgBD ) then !Hazelhos: added print statements
        write(*, '(a, 3i4)') dtxt//' LU loop i, j, iL: ', i, j, iL
        write(*,'(a,5es9.2)') dtxt//' Ra_ref, Rb, rsur, Ra_X, Ra_3m: ',&
                 L%Ra_ref, Rb, Rsur, L%Ra_X, L%Ra_3m
        write(*,'(a,6es9.2)') dtxt// &
          ' Vg_ref, Vg_eff, Gsto, Gigs, Gext:', Vg_ref,&  !DS skipped eff_fac here
            Vg_eff, Gsto, GigsBD, GextBD
     end if 

    !Hazelhos, autumn 2019: Xwater is added as input argument
    !QUERY Xemis does nothing! Why calculate it?

     call BiDirFluxes(L%is_water,L%SAI,Vg_ref,Gsto,&
                      GextBD,GigsBD,Xwater, &         !:inputs
                      Xtot,Xemis,debug=dbgBD)    !:outputs

     call CheckStop ( Xtot < 0.0, dtxt//'NEG XTOT'//LandDefs(iL)%name)


    ! Xtot is in ug(NH3)/m3.   V.X is in m/s * ug/m3 = ug/m2/s
    ! *14/17*3600 gives ugN/m2/h as wanted in Biogenics_mod
    ! QUERY - multiply by delta t?

    if ( dbgBD ) write(*,'(a, i3,4es12.3)') "BiDir_NHxemis 1", &
         iL, BiDir_2d%NHxEmis(i,j), Vg_ref, L%coverage, Xtot

    BiDir_2d%NHxEmis(i,j) = BiDir_2d%NHxEmis(i,j) + &
          Vg_ref * Xtot * L%coverage * (14./17*3600)
              !BUG? Vg_ref(icmp) * Xtot * L%coverage/360 * (14./17*3600)

    if ( dbgBD ) write(*,'(a,i3,4es12.3)') "BiDir_NHxemis 2",  &
        iL, BiDir_2d%NHxEmis(i,j), Vg_ref, L%coverage, Xtot

    BiDir_2d%Xtot(i,j) = BiDir_2d%Xtot(i,j) + Xtot * L%coverage

    ! Conc NH3 at 3-4m is needed:
    ! F = Vg(z)*C(z) = C(zmid)/(Ra+Rb+Rc) =  [C(zmid)-C(3m)]/Ra(zmid,3m)
    ! ie C(3m) = C(zmid) - F.Ra(zmid,3m)
    !BiDir
    ! F = Vg(z)*[C(zmid)-Cc] = [C(zmid)-C(3m)]/Ra(zmid,3m)
    ! ie C(3m) = C(zmid)-Vg.dRa*[C(zmid)-Cc] = C(zmid)[1-Vg.dRa] + Cc.Vg.dRa


    Ra_diff = L%Ra_ref - L%Ra_3m
    !Ra_diff = L%Ra_X - L%Ra_3m !TMP FIX QUERY!!!

    Ra_diff = max(Ra_diff, L%Ra_3m)
    if ( dbgBD ) write(*,'(a, 3es12.3)') "BiDir Ra_diff ", &
            L%Ra_ref, L%Ra_3m, Ra_Diff

    !FIXME if( Ra_diff < 0.0 ) Ra_diff = 0.0
    !FIXME  print '(a,9es10.3)', dtxt//'RADIFF', L%Ra_X, L%Ra_ref, L%Ra_3m, L%hveg, L%z0,Grid%z_ref,L%z_refd
    !call CheckStop(  Ra_diff < 0.0, dtxt// "BIDIRRa_diff")

    !JUL10 ugNH3_3m_LC   = ugNH3_zmid* (1- Ra_diff * Vg_ref )
    !JUL10 ugXH3_3m_LC   = Xtot *  Ra_diff * Vg_ref

    ugNH3_3m_LC   = ugNH3_zmid* (1- Ra_diff * Vg_eff ) ! Jul10 eff,

    !QUERY    Hazelhos 04-12-19: Why do we use Vg_eff here and
    ! not Vg_ref? at other places, which one we use depends on
    ! USES%EFFECTIVE_RESISTANCE (e.g. line 826)
    ! in Line 668, we use Ra_ref in stead of Ra_X. Ra_ref should be linked to Vg_ref?
    

    ugXH3_3m_LC   = Xtot *  Ra_diff * Vg_eff ! Jul10 eff

    if ( dbgBD ) write(*, '(a, i3,5es12.3)') "BiDir ug_NH3_grid 1 ", &
      iL, ugNH3_3m_grid, ugXH3_3m_grid, ugNH3_3m_LC, ugXH3_3m_LC, L%coverage

    !Hazelhos: For every land use, the concentrations are added depending on
    ! the coverage

    ugNH3_3m_grid = ugNH3_3m_grid + ugNH3_3m_LC * L%coverage
    ugXH3_3m_grid = ugXH3_3m_grid + ugXH3_3m_LC * L%coverage

    if ( dbgBD ) write(*, '(a, i3,3es12.3)') "BiDir ug_NH3_grid 2 ", &
      iL,  ugNH3_3m_grid, ugXH3_3m_grid, L%coverage

    if ( iL == 9 .and. Xtot > 45.0 ) then ! QUERY too high?
     write(*,'(a,2i3,2i4,es10.2,9f8.3)') dtxt//"XXBIG "//&
       trim(LandDefs(iL)%name), iL, current_date%hour, i_fdom(i), j_fdom(j), &
       GextBD, ugNH3_zmid, ugNH3_3m_LC, ugXH3_3m_LC, Xtot
    end if

    if( dbgBD .and. Xtot>25 ) then !  iL==1)  then
        call datewrite(dtxt//"Xtot>50:"//trim(LandDefs(iL)%name), [ L%SAI, &
          100*Gsto, oldRsur, Rsur, GextBD, ugNH3_zmid, ugNH3_3m_LC,Xtot, &
          ugXH3_3m_LC, BiDir_2d%NHxEmis(i,j) ],txt_pattern=.true.) !?? weird name
    end if

  end subroutine BiDir_ijFluxes
!----------------------------------------------------------------------------
  subroutine BiDir_ijFinish(i,j,gradfac,sumLand,DepLoss)
    !Last steps of Bi-Dir this i,j 
    !Sets gradient for NH3 and also SoilNH3 emission
    ! nb called from DryDep_module if Sumland > 0.01
    !QUERY Hazelhos 04-12-2019: What if Sumland = 0.02? The gridcell
    !  would be 98% sea but has no representation of it.
    !DS REPLY: we only use these gradient_facs to compare with observations,
    !which are assumed to be on land.

    integer, intent(in) :: i,j
    real, intent(inout) :: gradfac
    real, intent(in)    :: sumLand, DepLoss
    character(len=*), parameter :: dtxt='BiDirGradEm:'
    character(len=10) :: txt ! Land or Sea
    real :: ug_old

  ! 0) save values to BiDir fields

    ug_old             = BiDir_2d%NH3inst(i,j)
    BiDir_2d%NH3inst(i,j) = ugNH3_3m_grid     + ugXH3_3m_grid
    BiDir_2d%NH3_3m(i,j)  = BiDir_2d%NH3_3m(i,j) + ugNH3_3m_grid
    BiDir_2d%XH3_3m(i,j)  = BiDir_2d%XH3_3m(i,j) + ugXH3_3m_grid

    if(dbgBD) call datewrite(dtxt//"ugs: ", &
             -1, (/ ugNH3_zmid, ug_old, BiDir_2d%NH3inst(i,j) /) )

  ! 1) SoilNH3 emissions are used to set rcemis  ================

    SoilNH3(i,j) = BiDir_2d%NHxEmis(i,j)

    if ( dbgBD ) write(*,'(a, 3f7.3,3es12.3)') dtxt//" nh3: ", sumLand, ugNH3_zmid, gradfac, SoilNH3(i,j)

    if ( ugNH3_zmid < 1.0e-2 ) return ! Don't bother with tiny concs
    if ( sumLand <= 0.01     ) return ! No gradients over Sea. QUERY: RE-VISIT!

  ! 2) Sets new gradient of zmid/z3m =========================
    ! gradient_fac = ratio conc. at 3/50m. 

    gradfac = (ugNH3_3m_grid+ugXH3_3m_grid)/ugNH3_zmid

    if ( dbgBD ) then
       write(*,'(a,4(a,f10.4))') dtxt//'Land:', &
         ' ugNH3_3m', ugNH3_3m_grid, ' ugXH3_3m', ugXH3_3m_grid,&
         ' ugNH3_zmid', ugNH3_zmid,' gradfac', gradfac
    end if ! dbgBD
   
    ! TMP Just checking for extremes....
    if ( gradfac > 1.0 ) then

      if ( gradfac > 10.0 ) then
        write(*,'(a,2i4,5f9.3,es10.2)') dtxt//'XXGRAD '//trim(txt),&
          i_fdom(i),j_fdom(j), ugNH3_zmid, ugNH3_3m_grid, ugXH3_3m_grid,&
          BiDir_2d%Xtot(i,j), gradfac, BiDir_2d%NHxEmis(i,j)
      end if

     if (dbgBD) write(*,'(a,2f8.3,2es12.3)') dtxt//" NH3 Emis Dep Balance: ",&
         ugNH3_3m_grid, ugXH3_3m_grid,  BiDir_2d%NHxEmis(i,j), DepLoss

   end if
  end subroutine BiDir_ijFinish

!----------------------------------------------------------------------------
  subroutine BiDir_Derived(txt,n,limax,ljmax,nerr)
     character(len=*), intent(in) :: txt
     integer, intent(in) :: n,limax,ljmax
     integer, intent(inout) :: nerr ! set to -999 if problem
     integer :: i,j
     select case (txt)
       case("ugXTOT" )
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_INST) = BiDir_2d%Xtot(i,j)
         end forall
       case ( "ugNH3_3m" )
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_INST) = BiDir_2d%NH3_3m(i,j)
         end forall
      case ( "ugXH3_3m" )
        forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = BiDir_2d%XH3_3m(i,j)
        end forall
      !Not yet implemented. But Emis_mgm2_BioNatNH3 gives emissions
      !case ( "BiDir_NHx_Emissions" )
      !  forall ( i=1:limax, j=1:ljmax )
      !    d_2d( n, i,j,IOU_INST) = BiDir_2d%NHxEmis(i,j)
      !  end forall
      case default
        print *, "BiDir_Derived ERROR for"//trim(txt)
        nerr = -999
      end select

   end subroutine BiDir_Derived


!----------------------------------------------------------------------------
  ! quick function to replace YYYY MM etc with dates for load_input_data_BiDir
  function fixdates(fname) result (f)
    character(len=*), intent(in) :: fname
    character(len=200) :: f
     f=key2str(fname,'YYYY',iyr_trend)
     f=key2str(f,'MM',current_date%month )
     if ( current_date%month < 12 ) then 
       f=key2str(f,'ZZZZ',iyr_trend)
       f=key2str(f,'NN',current_date%month+1 )
     else
       f=key2str(f,'ZZZZ',iyr_trend+1)
       f=key2str(f,'NN', 1 ) !For Jan.
     end if
  end function  fixdates
!----------------------------------------------------------------------------
  subroutine load_input_data_BiDir()
     character(len=*), parameter :: dtxt='BIDIR:load monthly input data:'
     character(len=200) :: tmpdir !for SEASTUFF
     character(len=200) :: ifile  !for SEASTUFF
     integer :: i,j

     if ( current_date%month == old_month ) return
     if ( MasterProc) write(*,*) dtxt//'called DATE',step_main, iyr_trend, current_date 
 
     ! ============== MONTHLY UPDATES ========================
     !JUL2019 SEA STUFF
     ! 4) 
     ! NH4 is in mmole/m3 which is umole/L
     ! dates: use YYYY, MM and ZZZZ, NN (for end date if needed)

     ifile=trim(BiDir%InputDir)//'BiDirSea/FREEBIORYS2V4/YYYY/'// &
              'emep-ext-FREEBIORYS2V4_1mAV_YYYYMM01_ZZZZNN01_nh4.nc'
     ifile= fixdates(ifile)  ! replaces YYYY etc
       
     call ReadField_CDF(ifile, 'nh4', BiDir_2d%sea_nh4, nstart=1,&
       needed=.true.,debug_flag=.true.)
       !TMP interpol='zero_order',needed=.true.,UnDef=0.0, debug_flag=.true.)
       
     ifile=trim(BiDir%InputDir)//'BiDirSea/FREEBIORYS2V4/YYYY/'// &
              'emep-ext-FREEBIORYS2V4_1mAV_YYYYMM01_ZZZZNN01_ph.nc'
     ifile= fixdates(ifile)

     call ReadField_CDF(ifile, 'ph', BiDir_2d%sea_ph, nstart=1,&
       needed=.true.,debug_flag=.true.)

     if( debug_proc ) then
        write(*,*) 'TEST STR bidir',  trim(ifile), debug_li, debug_lj !
        write(*,*) 'TEST STR nh4', maxval(BiDir_2d%sea_nh4), minval(BiDir_2d%sea_nh4)
        write(*,*) 'TEST STR ph',  maxval(BiDir_2d%sea_ph), minval(BiDir_2d%sea_ph)
     end if

     ifile=trim(BiDir%InputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
        'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridT.nc'
     ifile= fixdates(ifile)

     call ReadField_CDF(ifile,'votemper', BiDir_2d%sea_gridT, nstart=1,&
       needed=.true.,debug_flag=.true.)

     !ifile=fix_dates( trim(BiDir%InputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
     !   'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridS.nc' )
     ifile= trim(BiDir%InputDir)//'BiDirSea/FREEGLORYS2V4/YYYY/'// &
        'emep-ext-GLORYS2V4_ORCA025_YYYYMM_gridS.nc' 
     ifile= fixdates(ifile)
     call ReadField_CDF(ifile,'vosaline', BiDir_2d%sea_gridS, nstart=1,&
       needed=.true.,debug_flag=.true.)

     if ( debug_proc ) write(*,*) minval(BiDir_2d%sea_nh4), minval(BiDir_2d%sea_gridT)
     call printCDF('BIDIR_SEA_nh4',BiDir_2d%sea_nh4,'seaNH4units')
     call printCDF('BIDIR_SEA_ph',BiDir_2d%sea_ph,'phUnits')
     call printCDF('BIDIR_SEA_gridT',BiDir_2d%sea_gridT,'seagridT')
     call printCDF('BIDIR_SEA_gridS',BiDir_2d%sea_gridS,'seagridS')

     old_month = current_date%month

  end subroutine load_input_data_BiDir
  
  
  subroutine set_BiDirSea(i,j,ssNH4,sspH,sstK,S)
    integer, intent(in) :: i,j
    real, intent(out) :: sstK, ssNH4, sspH, S
     ssNH4  = BiDir_2d%sea_nh4(i,j)
     sspH   = BiDir_2d%sea_ph(i,j)
     sstK = BiDir_2d%sea_gridT(i,j) + 273.15  ! or sstK_nwp  TO DO
     S    = BiDir_2d%sea_gridS(i,j)
  end subroutine set_BiDirSea

end module BiDir_emep
