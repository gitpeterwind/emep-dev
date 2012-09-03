! < ForestFire_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
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

module ForestFire_ml
 !----------------------------------------------------------------
 ! Uses emissions from either:
 ! 1) FINNv1 daily data 2002 - 2011
 ! REFERENCE:
 ! Wiedinmyer, C., Akagi, S. K., Yokelson, R. J., Emmons, L. K., Al-Saadi,
 ! J. A., Orlando, J. J., and Soja, A. J.: The Fire INventory from NCAR (FINN) 
 ! - a high resolution global model to estimate the emissions from open 
 !  burning, Geosci. Model Dev. Discuss., 3, 2439-2476, 
 !   doi:10.5194/gmdd-3-2439-2010, 2010.
 ! http://www.geosci-model-dev-discuss.net/3/2439/2010/gmdd-3-2439-2010.html
 !
 ! 2)  GFED 3 (Global Forest Emission database)
 ! http://www.falw.vu/~gwerf/GFED/
 ! Currently programmed for 8-daily data (available for 2001 - 2007)
 !
 ! 3) GFASv1 Real-Time Fire Emissions
 ! Daily data. Available since 2003 or 2011, depending on version, from MARS
 !   http://www.gmes-atmosphere.eu/about/project_structure/input_data/d_fire/ProductsInMARS/
 ! REFERENCE:
 ! Kaiser, J.W., Heil, A., Andreae, M.O., Benedetti, A., Chubarova, N.,
 !   Jones, L., Morcrette, J.-J., Razinger, M., Schultz, M. G., Suttie, M.,
 !   and van der Werf, G. R.: Biomass burning emissions estimated with a global
 !   fire assimilation system based on observed fire radiative power,
 !   Biogeosciences, 9, 527-554, doi:10.5194/bg-9-527-2012, 2012.
 !----------------------------------------------------------------
  use CheckStop_ml,      only : CheckStop
  use ChemChemicals_ml,  only : species
  use ChemSpecs_tot_ml
  use EmisDef_ml,        only : ISNAP_NAT &! Fires assigned to SNAP-11 usually
   ,NEMIS_FILE, EMIS_FILE  ! which pollutants are wanted, e.g. sox, pm25

  use GridValues_ml,     only : i_fdom, j_fdom, debug_li, debug_lj, &
                                 debug_proc,xm2,GRIDWIDTH_M
  use Io_ml,             only : PrintLog, datewrite
  use MetFields_ml,      only : z_bnd
  use ModelConstants_ml, only : MasterProc, KMAX_MID, &
                                USE_FOREST_FIRES, DEBUG_FORESTFIRE, FORECAST, &
                                IOU_INST,IOU_HOUR,IOU_HOUR_MEAN, IOU_YEAR
  use NetCDF_ml,         only : ReadField_CDF, Out_netCDF,Real4 ! Reads, writes 
  use OwnDataTypes_ml,   only : Deriv, TXTLEN_SHORT
  use Par_ml,            only : MAXLIMAX, MAXLJMAX, &
                                  me,limax,ljmax
  use PhysicalConstants_ml, only : AVOG
  use Setup_1dfields_ml, only : rcemis
  use SmallUtils_ml,     only : find_index
 ! No. days per year, date-type :
  use TimeDate_ml,only : nydays, nmdays, date, current_date  
  use TimeDate_ExtraUtil_ml,  only: date2string
  implicit none

!  Unimod calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!

  public :: Fire_Emis
  public :: Fire_rcemis
  private :: Export_FireNc

  logical, public, allocatable, dimension(:,:), save ::  burning
  real, private, allocatable, dimension(:,:,:), save :: BiomassBurningEmis

  integer, private, save ::  ieCO=-1 ! index for CO

  character(len=TXTLEN_SHORT), private :: FF_poll
  integer :: iemep

  !/ Defintions of GFED data. If known, we assign the GFED pollutant which
  !  corresponds to each possible EMEP emission file. Simply add EMEP
  !  lines as required - be consistent with EmisDefs though. (We can
  !  have more definitions than used in EmisDefs, but not vice.versa.

  ! Assign mol. wts of the GFED data  where known. If mol. wt set to
  ! zero, the code in Fire-rcemis will use the values from the
  ! ChemSpecs_ml, species()%molwt.
  !
  ! If GFED doesn't have emissionss, set a "-" for GFED, then the
  ! desired emission factor (g/kg DW), and then follow the
  ! example in Fire_setups for NH3:
  !

  ! =======================================================================
  !  Mapping to EMEP species
  !
  type, private :: bbtype
    character(len=TXTLEN_SHORT) :: txt
    real :: unitsfac
    real :: frac
    integer :: emep
  endtype bbtype

  include 'BiomassBurningMapping.inc'

  ! just to keep track
  real, private, save ::  sum_emis(NBB_DEFS) ! debug
  ! matrix to get from forest-fire species to EMEP ones

  type, private :: bbmap
    integer :: emep
    integer :: bbspec
    real    :: fac
  end type bbmap
  type(bbmap), private, dimension(NCMSPECS) :: fmap

  ! =======================================================================


contains
subroutine Fire_Emis(daynumber)
!.....................................................................
!**    DESCRIPTION:
!    Reads forest-fire emissions. So far set up for GFED 8d, but in
!    principal we can re-code by simply adding alternative
!    subroutines, e.g. to cope with pre-2001 monthly emissions

  integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


  real,allocatable :: rdemis(:,:)  ! Emissions read from file
  integer :: i,j,nstart, alloc_err
  logical, save :: my_first_call = .true.   ! DSFF
  logical :: my_first_defs = .true. 
  integer :: dd_old = -1,  n
  real    :: fac, to_kgm2s   
  character(len=TXTLEN_SHORT), dimension(NBBSPECS) :: FF_names

  logical, save :: calc_remppm = .false.      ! Needed to get REMPPM25
  integer, save :: ieOC=-1,ieBC=-1,iePM25=-1  !    " / "
  real    :: OMbb
  integer :: loc_maxemis(2) ! debug

  character(len=*), parameter :: &
    GFED_PATTERN = 'GFED_ForestFireEmis.nc',&
    FINN_PATTERN = 'FINN_ForestFireEmis_YYYY.nc',&
    GFAS_PATTERN = 'GFAS_ForestFireEmis_YYYY.nc'
  character(len=len(GFAS_PATTERN)) :: fname = ''
  logical :: my_debug=.false.
  integer, parameter :: verbose = 1

  if(my_first_call) &
    call PrintLog("Biomass Mapping: "//trim(BiomassBurningMapping),MasterProc)

  select case(verbose)
    case(:0);my_debug=.false.
    case(1) ;my_debug=DEBUG_FORESTFIRE.and.MasterProc.and.my_first_call
    case(2) ;my_debug=DEBUG_FORESTFIRE.and.MasterProc
    case(3) ;my_debug=DEBUG_FORESTFIRE
    case(4:);my_debug=.true.
  endselect

  nstart = -1 ! reset for GFED
  select case(BiomassBurningMapping(1:4))
  case("GFED") ! 8-day values
    select case(current_date%year)
    case(2001:2007)
      if(MasterProc)&
        write(*,*) "WARNING! FFIRE GFED USED! May not be working properly check results!"
    case default
      if(my_first_call)&
        call PrintLog("8d GFED Forest Fires: only between 2001--2007",MasterProc)
      my_first_call = .false.
      return
    endselect
    if(DEBUG_FORESTFIRE.and.MasterProc) &
      write(*,*) "Into the FIRE days:", current_date%year, &
        daynumber, dd_old, mod(daynumber,8), my_first_call

    ! GFED Fire emissions are called at 8 days intervals (1, 9, 17, ....)
    ! 46 values available each year: day 361 is the last one.
    ! Return unless new period

    if(.not.my_first_call.and.mod(daynumber,8)/= 1) return
    nstart=(current_date%year-2001)*46+(daynumber+7)/8
  case("FINN")
  case("GFAS")
  case default
    call CheckStop("Unknown B.B.Mapping: "//trim(BiomassBurningMapping))
  endselect

  if(DEBUG_FORESTFIRE.and.MasterProc) &
    write(*,*) "Into the FIRE days:", current_date%year, &
      daynumber, dd_old, mod(daynumber,8), my_first_call

  if(dd_old==daynumber) return   ! Only calculate once per day max
  dd_old = daynumber

  if(my_first_call)then
    call Fire_setup()   ! Gets InvMolwWtFac

    ! Assume NEMIS_FILE for now
    allocate(BiomassBurningEmis(NBBSPECS,MAXLIMAX,MAXLJMAX),&
             burning(MAXLIMAX,MAXLJMAX),stat=alloc_err)
    call CheckStop(alloc_err,"ForestFire BiomassBurningEmis alloc problem")
    BiomassBurningEmis(:,:,:)=0.0
    my_first_call = .false.

    do n=1,NBBSPECS
      select case(species(FF_defs(n)%emep)%name)
      ! CO is special. Keep the index
      case("CO"            );ieCO  =n
      ! GFED & FINN have OC and BC as well as PM25.
      ! Where the emep model asks for FFIRE_REMPPM25 we need to get this by substraction
      case("FFIRE_REMPPM25");iePM25=n;calc_remppm=.true.
      case("FFIRE_OC"      );ieOC  =n
      case("FFIRE_BC"      );ieBC  =n
      endselect
    enddo
    call CheckStop(ieCO<0,"No mapping for 'OC' found on "//BiomassBurningMapping)
    if(calc_remppm)then
      call CheckStop(ieOC<0,"No mapping for 'FFIRE_OC' found on "//BiomassBurningMapping)
      call CheckStop(ieBC<0,"No mapping for 'FFIRE_BC' found on "//BiomassBurningMapping)
    endif
  endif
  allocate(rdemis(MAXLIMAX,MAXLJMAX),stat=alloc_err)
  call CheckStop(alloc_err,"ForestFire rdemis alloc problem")

  if(DEBUG_FORESTFIRE .and. MasterProc) &
    write(*,*) "FOREST_FIRE: ", daynumber,nstart, NBBSPECS

  ! We need to look for forest-fire emissions which are equivalent
  ! to the standard emission files:

  do n = 1, NBBSPECS
    FF_poll = trim(FF_defs(n)%txt)
    iemep = FF_defs(n)%emep

    if(DEBUG_FORESTFIRE.and.MasterProc) &
      write(*,"( a,2i3,1x,2a8,2i3,a)") "FIRE SETUP: ", &
        n, iemep, trim(species(iemep)%name), trim(FF_poll), len_trim(FF_poll)
! FORECAST mode: if file/variable/timestep not found it should not crash
    rdemis(:,:)=0.0


    select case(BiomassBurningMapping(1:4))
    case("GFED")
      fname = date2string(GFED_PATTERN,current_date)
      if(my_debug) &
        write(*,*) "FFIRE GFED ", me, n, nstart,  trim(FF_poll), trim(fname)
      call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol='zero_order',&
        needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG_FORESTFIRE)
      !unit conversion to GFED [g/m2/8day]->[kg/m2/s]
      to_kgm2s = 1.0e-3 /(8*24.0*3600.0)
!     rdemis = rdemis * to_kgm2s
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*to_kgm2s
    case("FINN")
      fname = date2string(FINN_PATTERN,current_date)
      if(my_debug) &
        write(*,*) "FFIRE FINN ", me, n, daynumber,  trim(FF_poll), trim(fname)
      call ReadField_CDF(fname,FF_poll,rdemis,daynumber,interpol='mass_conservative',&
        needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG_FORESTFIRE)
      !unit conversion to FINN:
      fac=FF_defs(n)%unitsfac                       ! --> [kg or mole /day]
      fac=fac/(GRIDWIDTH_M*GRIDWIDTH_M*24.0*3600.0) ! [kg or mole /day]->[kg/m2/s]
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac*xm2(i,j)
    case("GFAS")
      fname = date2string(GFAS_PATTERN,current_date)
      nstart = daynumber
! something more sophisticated is needed for YYYY_ or YYYYMM_ files,
! e.g. use ReadTimeCDF and nctime2idate/idate2nctime to find the right record:
!  nstart=FindTimeCDFRecord(fname,current_date,prec_ss=3600.0*12)
      if(my_debug) &
        write(*,*) "FFIRE GFAS ", me, n, nstart,  trim(FF_poll), trim(fname)
      call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol='conservative',&
          needed=.not.FORECAST,debug_flag=DEBUG_FORESTFIRE)
      ! GFAS units are [kg/m2/s]. No unit conversion is needed.
    endselect

    ! Assign . units should be [kg/m2/s] here
   !BiomassBurningEmis(n,:,:) = rdemis(:,:)
    forall(j=1:ljmax,i=1:limax) BiomassBurningEmis(n,i,j) = rdemis(i,j)

    if(my_first_defs) call PrintLog(&
      "ForestFire_ml :: Assigns "//trim(FF_poll) , MasterProc)

    if(DEBUG_FORESTFIRE) sum_emis(n)=sum_emis(n)+sum(BiomassBurningEmis(n,:,:))
  enddo ! BB_DEFS
  my_first_defs  = .false.
  deallocate(rdemis)

  ! Logical to let Unimod know if there is any emission here to worry about
  burning(:,:) = ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )

  ! If needed, calculate REMPPM25, assuming OM/OC = 1.3  
  !RB: But we assume OM/OC=1.7 for forest fire emissions in the EMEP model? Test changing to 1.7 here also
  if(calc_remppm) then
    do j=1,ljmax
      do i=1,limax
!AMVB:  OMbb = 1.7*BiomassBurningEmis(ieOC,i,j) + BiomassBurningEmis(ieBC,i,j)
!AMVB:  This is already taken into account by FF_defs(ieOC)%unitsfac      
!BM     OMbb=BiomassBurningEmis(ieOC,i,j)+BiomassBurningEmis(ieBC,i,j)
        OMbb=1.7*BiomassBurningEmis(ieOC,i,j)+BiomassBurningEmis(ieBC,i,j)
        if(DEBUG_FORESTFIRE.and.MasterProc.and.burning(i,j)) &
          write(*,"(a,2i4,3es10.3)") "BURN REMPPM25, ", i_fdom(i), j_fdom(j), &
            BiomassBurningEmis(iePM25,i,j),OMbb,BiomassBurningEmis(iePM25,i,j)-OMbb
        BiomassBurningEmis(iePM25,i,j) = DIM(BiomassBurningEmis(iePM25,i,j),OMbb)
          !max(BiomassBurningEmis(iePM25,i,j)-OMbb,0.0)
      enddo
    enddo
  endif

  if(DEBUG_FORESTFIRE.and.debug_proc) then
    n = ieCO
    FF_poll = trim(FF_defs(n)%txt)
    loc_maxemis = maxloc(BiomassBurningEmis(n,:,: ) )
    call datewrite("SUM_FF: "//trim(FF_poll),  &
      (/ daynumber, n, i_fdom(loc_maxemis(1)), j_fdom(loc_maxemis(2)) /) ,&
      (/  sum_emis(n), maxval(BiomassBurningEmis(n,:,: ) ), &
          BiomassBurningEmis(n,debug_li,debug_lj) /) &
    )
  endif ! debug_proc
endsubroutine Fire_Emis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Fire_setup()
! Pre-calculate conversion factors to get
! from BiomassBurning's kg/m2/s to molecules/cm3/s.

  integer :: idef, f, n, alloc_err, iFF, nFF
  integer :: iemep, nemep
  character(len=TXTLEN_SHORT) :: last_FF_poll = "NotSet"
  real :: fac

  nFF   = 0
  nemep = 0

  do n = 1, NBB_DEFS
    iemep   = FF_defs(n)%emep
    nemep   = nemep  + 1
    FF_poll = trim(FF_defs(n)%txt)
    if ( FF_poll /= last_FF_poll ) then
      nFF = nFF + 1
      last_FF_poll = trim(FF_poll)
!     if(MasterProc) print *, "MATCH+1 FFIRE ", n, iemep, nemep, nFF, trim(FF_defs(n)%txt)
    endif
!   if(MasterProc) print *, "END MATCH FFIRE ", n, iemep, nemep, nFF, trim(FF_defs(n)%txt)

  !// And one final conversion factor.
  ! The biomassBurning array is kept in kg/m2/s for consistency with other
  ! emissions. We will later convert to molecules/cm3/s after spreading
  ! through a vertical distance dz
  !
  ! If we had E in kg/m2/s, we would then take
  !  E*1.0e3  -> g/m2/s
  !  E*0.1    -> g/cm2/s
  !  E*0.1 /MW * Av -> molec/cm2/s
  !  E*0.001 /MW * Av / DZ -> molec/cm3/s where DZ is spread in m
  !  i.e. fmap should be 0.001*Av/MW
  !  (plus account for the fraction of the inventory assigned to EMEP species)

    fac =  &!BUG  FF_defs(n)%unitsfac  &  ! MW scale if needed
          FF_defs(n)%frac      &  ! mass-fraction split
        * 0.001 * AVOG /species(iemep)%molwt

    fmap(nemep) = bbmap( iemep,  nFF, fac ) 
!       if ( MasterProc ) print *, "FFIRE MAP ", nemep, fmap(nemep)
  enddo !  n
endsubroutine Fire_setup
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Fire_rcemis(i,j)
!  Disperses the fire emissions vertically and converts to molecules/cm3/s.

!// Injection height: here over 8 levels. Alternative could be PBL
!   or  equally upto ca. 2*PBL (suggested by Sofiev, GEMS)
!ds QUERY - should the emissions be divided equally by level?
!   - will give a higher mixing ratio for thinner levels

  integer, intent(in) :: i,j

  integer, parameter :: KEMISFIRE = 12
  real, dimension(KEMISFIRE:KMAX_MID) :: invDeltaZfac !  height of layer in m div 9
  integer ::  k, n, iem, iFF

  integer, parameter ::  N_LEVELS = KMAX_MID - KEMISFIRE + 1  ! = 9.0 here

  real    :: origrc, bbe
  logical :: debug_flag

  debug_flag = (DEBUG_FORESTFIRE.and.debug_proc .and.&
                i==debug_li.and.j==debug_lj)
  if(debug_flag.and.BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
    write(*,"(a,5i4,es12.3,f9.3)") "BurningDEBUG ", me, i,j, &
      i_fdom(i), j_fdom(j), BiomassBurningEmis(ieCO,i,j)

  !/ Here we just divide by the number of levels. Biased towards
  !  different levels since thickness and air content differ. Simple though.

  do k = KEMISFIRE, KMAX_MID
    invDeltaZfac(k) = 1.0/ (z_bnd(i,j,k) - z_bnd(i,j,k+1)) /N_LEVELS
  enddo
 
 do n = 1, NCMSPECS
    iFF = fmap(n)%bbspec
    iem = fmap(n)%emep
    bbe = BiomassBurningEmis(iFF,i,j) * fmap(n)%fac
    origrc = rcemis(iem,KMAX_MID)   ! just for printout

    ! distribute vertically:
    do k = KEMISFIRE, KMAX_MID
      rcemis(iem,k) = rcemis(iem,k) + bbe*invDeltaZfac(k)
    enddo !k

    if(debug_flag) then
      k=KMAX_MID
      write(*,"(a,2i3,1x,a8,i4,es10.2,4es10.2)") "FIRERC ",&
        n, IFF, trim(species(iem)%name), k, BiomassBurningEmis(iFF,i,j),&
        fmap(n)%fac, invDeltaZfac(k), origrc, rcemis( iem, k )
    endif

!DSBB    !--  Add up emissions in ktonne ......
!DSBB    !   totemadd(iem) = totemadd(iem) + &
!DSBB    !   tmpemis(iqrc) * dtgrid * xmd(i,j)

  enddo ! n
 !       call Export_FireNc() ! Caused problems on last attempt

endsubroutine Fire_rcemis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Export_FireNc()
  type(Deriv) :: def1 ! definition of fields
 
  def1%class='ForestFireEmis' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0        !not used
  def1%name='CO'        !written
  def1%unit='g/m2'      !written
  def1%name='CO_zero'
  def1%name='CO_ASCII'       

  call Out_netCDF(IOU_INST,def1,2,1, BiomassBurningEmis(ieCO,:,:),1.0,&
                  CDFtype=Real4,fileName_given='FF.nc')
endsubroutine Export_FireNc

endmodule ForestFire_ml
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
