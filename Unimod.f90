! <Unimod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007 met.no
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
program myeul
  !
  !     This is the main program for the off-line regional scale multilayer
  !     eulerian model at emep/msc-w. the main program contains the outer
  !     time-loop which runs through all time-levels a new meteorological
  !     data-set is read into the model from file. the inner time-loop
  !     runs through the physical time-step.
  !
  !------------------------------------------------------------------------

  use My_Outputs_ml,    only : set_output_defs &
                              ,NHOURLY_OUT
  use My_Timing_ml,     only : lastptim,mytimm,Output_timing, &
       Init_timing, Add_2timing, Code_timer, &
       tim_before,tim_before0,tim_before1, &
       tim_after,tim_after0

  use Advection_ml,     only : vgrid,adv_var, assign_nmax,assign_dtadvec
  use Aqueous_ml,       only : init_aqueous, Init_WetDep   !  Initialises & tabulates
  use AirEmis_ml,       only : aircraft_nox, lightning
  use Biogenics_ml,     only : Init_BVOC, SetDailyBVOC  !dsBVOC
! hb NH3emis
  use calc_emis_potential_ml,only : NH3emis_potential,lNH3emis_pot, &
       readNH3emis,lEmis50_nh3
  use NH3Emis_variation_ml,  only : SetNH3 
  use BoundaryConditions_ml, only : BoundaryConditions
  use CheckStop_ml,     only : CheckStop
  use ChemChemicals_ml, only : define_chemicals
  use ChemGroups_ml,    only : Init_ChemGroups  
  use DefPhotolysis_ml, only : readdiss
  use Derived_ml,       only : Init_Derived
  use DerivedFields_ml, only : f_2d, f_3d
  use DO3SE_ml,         only : Init_DO3SE !LPJ
  use EcoSystem_ml,     only : Init_EcoSystems
  use EmisDef_ml,       only : AIRNOX, NH3EMIS_VAR ! NH3emis, experimental
  use Emissions_ml,     only : Emissions ,newmonth      !  subroutines
  use ForestFire_ml,    only : Fire_Emis
  use GridValues_ml,    only : MIN_ADVGRIDS,GRIDWIDTH_M,Poles
  use Io_ml  ,          only : IO_MYTIM,IO_RES,IO_LOG,IO_TMP,IO_DO3SE,&
                               IO_NH3_DEB ! hb NH3Emis
  use Io_Progs_ml  ,    only : read_line, PrintLog
  use Landuse_ml,       only : InitLandUse, SetLanduse, Land_codes
  use MassBudget_ml,    only : Init_massbudget,massbudget
  use Met_ml,           only : metvar,MetModel_LandUse,&
       Meteoread,MeteoGridRead
  use ModelConstants_ml,only : KMAX_MID   &   ! No. vertical layers
       ,MasterProc &   ! set true for host processor, me==0
       ,RUNDOMAIN  &   !
       ,NPROC      &   ! No. processors
       ,METSTEP    &   ! Hours between met input
       ,runlabel1  &   ! explanatory text
       ,runlabel2  &   ! explanatory text
       ,nprint,nterm,iyr_trend, PT &
       ,IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY &
       ,USE_CONVECTION, USE_SOILWATER,USE_SOIL_NOX,USE_FOREST_FIRES &
       ,USE_BVOC_2010,USE_DUST,DO_SAHARA &
       ,NBVOC  & !
       ,FORECAST       ! FORECAST mode
  use NetCDF_ml,        only : Init_new_netCDF
  use OutputChem_ml,    only : WrtChem
  use Par_ml,           only : me,GIMAX,GJMAX ,MSG_MAIN1,MSG_MAIN2&
       ,Topology, parinit &
       ,li0,li1,lj0,lj1,tgi0,tgj0  &   ! FOR TESTING
       ,limax, ljmax, MAXLIMAX, MAXLJMAX, gi0, gj0
  use PhyChem_ml,       only : phyche  ! Calls phys/chem routines each dt_advec
  use Sites_ml,         only : sitesdef  ! to get output sites
  use Tabulations_ml,   only : tabulate
  use TimeDate_ml,      only : date, current_date, day_of_year,daynumber,&
                               startdate,enddate
  use TimeDate_ExtraUtil_ml,only : assign_NTERM,date2string
  use Trajectory_ml,    only : trajectory_init,trajectory_in
  use Nest_ml,          only : wrtxn &
       ,FORECAST_NDUMP &  ! number of nested outputs on FORECAST mode
       ,outdate           ! dates of nested outputs on FORECAST mode
  !--------------------------------------------------------------------
  !
  !  Variables. There are too many to list here. Still, here are a
  !  few key variables that  might help:

  !     dt_advec       - length of advection (phyche) time-step
  !     GRIDWIDTH_M    - grid-distance
  !     gb             - latitude (sorry, still Norwegian influenced..)
  !     NPROC          - number of processors used
  !     me             - number of local processor, me=0 is host (=MasterProc)
  !                       processor where many read/writes are done
  !     ndays          - number of days since 1 january (max 365 or 366)
  !     thour          - utc-time in hours every time-step
  !     .
  !     .
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !
  !
  !
  !     +---------------------------------------------------------+
  !     +                                                         +
  !     +                                                         +
  !     +                  start main programme                   +
  !     +                                                         +
  !     +_________________________________________________________+
  !
  !
  !
  !     declarations in main programme.
  !
  implicit none

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO

  logical, parameter :: DEBUG_UNI = .false.
  integer numt, nadd, oldseason,newseason
  integer i
  integer :: mm, mm_old   ! month and old-month
  integer :: nproc_mpi,cyclicgrid
  character (len=230) :: errmsg,txt

  !
  !     initialize the parallel topology
  !
  nproc_mpi = NPROC
  CALL MPI_INIT(INFO)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_mpi, INFO)

  !  some checks
  !
  if(nproc_mpi /= NPROC)then
     if(me==0)write(*,fmt=*)"Wrong processor number!", &
          " Program was compiled with NPROC = ",NPROC, &
          " but MPI found ", nproc_mpi," processors available.", &
          " Please change NPROCX or NPROCY in ModelConstants_ml.f90"
          CALL MPI_FINALIZE(INFO)
          stop
  end if
  call CheckStop( digits(1.0) < 50, &
       "COMPILED WRONGLY: Need double precision, e.g. f90 -r8")

  call parinit(MIN_ADVGRIDS)     !define MasterProc subdomains sizes and position

  if (MasterProc) then
     open(IO_RES,file='eulmod.res')
     open(IO_LOG,file='RunLog.out')
     open(IO_TMP,file='INPUT.PARA')
  endif

  call read_line(IO_TMP,txt,status(1))
  read(txt,*) iyr_trend

  call read_line(IO_TMP,runlabel1,status(1))! explanation text short
  call read_line(IO_TMP,runlabel2,status(1))! explanation text long
  call read_line(IO_TMP,txt,status(1))  ! meteo year,month,day to start the run
  read(txt,*)startdate(1:3)             ! meteo hour to start the run is set in assign_NTERM
  call read_line(IO_TMP,txt,status(1))  ! meteo year,month,day to end the run
  read(txt,*)enddate(1:3)               ! meteo hour to end the run is set in assign_NTERM

  if(FORECAST)then  ! read dates of nested outputs on FORECAST mode
    do i=1,FORECAST_NDUMP
      call read_line(IO_TMP,txt,status(1))
      read(txt,"(I4,3(I2),I4)")outdate(i)  ! outputdate YYYYMMDDhhssss
      if(MasterProc) print "(1X,A)",&
        date2string("Forecast nest/dump at YYYY-MM-DD hh:mm:ss",outdate(i))
    enddo
  end if

  if( MasterProc ) then
    close(IO_TMP)
    call PrintLog( trim(runlabel1) )
    call PrintLog( trim(runlabel2) )
    call PrintLog( date2string("startdate = YYYYMMDD",startdate(1:3)) )
    call PrintLog( date2string("enddate   = YYYYMMDD",enddate  (1:3)) )
    write(unit=txt,fmt="(a,i4)") "iyr_trend= ", iyr_trend
    call PrintLog( trim(txt) )
    write(unit=IO_LOG,fmt="(a12,4i4)")"RunDomain:  ", RUNDOMAIN

   ! And record some settings to RunLog (will recode later)
    if(  FORECAST       ) call PrintLog("Forecast mode on")
    call PrintLog("Options used of (convec., soilwater, soilnox, forest fires)")
    if(  USE_CONVECTION ) call PrintLog("Convection used")
    if(  USE_SOILWATER  ) call PrintLog("SoilWater  switch on")
    if(  USE_SOIL_NOX   ) call PrintLog("SoilNOx    switch on")
    if(  USE_FOREST_FIRES)call PrintLog("ForestFires switch on")
    if(  USE_BVOC_2010 )  call PrintLog("BVOC 2010 switch on")
    call PrintLog("Options used of (dust, sahara)")
    if(  USE_DUST        )call PrintLog("Dust switch on")
    if(  DO_SAHARA       )call PrintLog("Sahara switch on")
endif



  !*** Timing ********

  call Init_timing()
  call Code_Timer(tim_before0)
  tim_before = tim_before0

  call MeteoGridRead(cyclicgrid) !define grid projection and parameters
  call Topology(cyclicgrid,Poles)!def GlobalBoundaries & subdomain neighbors
  call assign_NTERM(NTERM)!set NTERM, the number of 3-hourly periods
  call assign_dtadvec(GRIDWIDTH_M)! set dt_advec



  !     Decide the frequency of print-out
  !
  nadd = 0
  nprint = nterm
  if (nterm >  nprint) nadd = 1

  if (me ==  0) write(6,*)'nterm, nprint',nterm, nprint

  !-------------------------------------------------------------------
  !
  !++  parameters and initial fields.
  !

  call Add_2timing(1,tim_after,tim_before,"Before define_Chemicals")

  call define_chemicals()    ! sets up species details
  call Init_ChemGroups()    ! sets up species details

  call assign_nmax(METSTEP)   ! No. timesteps in inner loop

  call trajectory_init()

  call Add_2timing(2,tim_after,tim_before,"After define_Chems, readpar")


  call MeteoRead(1)

  call Add_2timing(3,tim_after,tim_before,"After infield")

  if ( MasterProc .and. DEBUG_UNI) write(6,*)"Calling emissions with year" ,&
          current_date%year

  call Emissions(current_date%year)

! hb NH3Emis
!    if (NH3EMIS_VAR) then                                                     
! writes temporal emis variation for Tange to file                             
!       open(IO_NH3_DEB,FILE='out.Tange.dat')                                  
!       write(IO_NH3_DEB,'(4a7,18a12)')"mm","dd","hh","TIME1","ISO_STABLE",&   
!                      "OPEN_STABLE","STORAGE","WIN_CROP","SPR_CROP",&       
!                      "SPR_SBEET","SPR_GRASS","MANURE1","MANURE2","MANURE3",& 
!                      "MANURE4","MANURE4a","MIN_SPRING","MIN_AUTUMN",&      
!                      "GRAZ_CATTLE","NH3_GRASS","TRAFFIC","SUM "            
!    endif       
  if (NH3EMIS_VAR)call readNH3emis() !read 16.7km activity NH3 emissions     
  if (NH3EMIS_VAR)write(6,*)'New ammonia emissions on proc ', &         
        me,sum( lEmis50_nh3(:,:,:) )

  if (NH3EMIS_VAR)call NH3emis_potential(current_date%year) !calc emission potential                                                                      
                                                                               
  if (NH3EMIS_VAR)write(6,*)'Potential emissions on proc ', &
          me,sum( lNH3emis_pot(:,:,:) )


  ! daynumber needed  for BCs, so call here to get initial

  daynumber=day_of_year(current_date%year,current_date%month,current_date%day)


  call MetModel_LandUse(1)   !

  call InitLandUse()  !  Reads Inputs.Landuse, Inputs.LandPhen

! Read data for DO3SE (deposition O3 and  stomatal exchange) module
! (also used for other gases!)

        call Init_DO3SE(IO_DO3SE,"Inputs_DO3SE.csv",Land_codes, errmsg)
        call CheckStop(errmsg, "Reading DO3SE ")

  call Init_EcoSystems()     !   Defines ecosystem-groups for dep output

  call Init_Derived()        ! Derived field defs.

  call Init_BVOC()

  call tabulate()    ! =>  sets up tab_esat, etc.

  call Init_WetDep()   ! sets up scavenging ratios

  call set_output_defs()     ! Initialises outputs
  call sitesdef()      !--- see if any output for specific sites is wanted
                       !   (read input files "sites.dat" and "sondes.dat" )

  call vgrid    !  initialisation of constants used in vertical advection
  if ( MasterProc .and. DEBUG_UNI) write(6,*)"vgrid finish"

! open only output netCDF files if needed
  if ( MasterProc ) then

    print *, "NETCDFINITS; maxval", maxval(f_2d%iotype)
    call Init_new_netCDF(trim(runlabel1)//'_inst.nc',IOU_INST)
    call Init_new_netCDF(trim(runlabel1)//'_fullrun.nc',IOU_YEAR)

    if( maxval(f_2d%iotype) == IOU_HOUR ) &
      call Init_new_netCDF(trim(runlabel1)//'_hour.nc',IOU_HOUR)
    if( maxval(f_2d%iotype) >= IOU_DAY ) &
      call Init_new_netCDF(trim(runlabel1)//'_day.nc',IOU_DAY)
    if( maxval(f_2d%iotype) >= IOU_MON ) &
      call Init_new_netCDF(trim(runlabel1)//'_month.nc',IOU_MON)

!FEB2011 
    !if(any(f_2d%inst).or.any(f_3d%inst))&
    !  call Init_new_netCDF(trim(runlabel1)//'_inst.nc',IOU_INST)
    !!netCDF hourly is initiated in Output_hourly
!
!FEB2011 - keep with NHOURLY test for now
    if(NHOURLY_OUT>0)&
      call Init_new_netCDF(trim(runlabel1)//'_hour.nc',IOU_HOUR)
!
!    if(any(f_2d%day).or.any(f_3d%day))&
!      call Init_new_netCDF(trim(runlabel1)//'_day.nc',IOU_DAY)
!
!    if(any(f_2d%month).or.any(f_3d%month))&
!      call Init_new_netCDF(trim(runlabel1)//'_month.nc',IOU_MON)

    ! The fullrun file contains the accumulated or average results
    ! over the full run period, often a year, but even just for
    ! a few timesteps if that is all that is run:
    !FEB2011 if(any(f_2d%year).or.any(f_3d%year))&
    !FEB2011 - Always....  call Init_new_netCDF(trim(runlabel1)//'_fullrun.nc',IOU_YEAR)
  endif

  call metvar(1)

  call adv_var(1)

  call Add_2timing(4,tim_after,tim_before,"After tabs, defs, adv_var")
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !
  !     performance of physical and chemical calculations
  !
  !     three-hourly time loop starts here.
  !

  mm_old = 0
  oldseason = 0

  tim_before = tim_before0
  call Add_2timing(5,tim_after,tim_before,"Total until numt loop")
  call Code_timer(tim_before1)


  do numt = 2, nterm + nadd

! hb NH3emis                                                                
     if (NH3EMIS_VAR) call SetNH3()

     mm = current_date%month

     if (mm == 12 .or. mm < 3) then
        newseason = 1
     else if(mm < 6) then
        newseason = 2
     else if(mm < 9) then
        newseason = 3
     else if(mm < 12) then
        newseason = 4
     end if

     ! - daynumber needed  for BCs, so call here to be safe

     daynumber=day_of_year(current_date%year,current_date%month,&
          current_date%day)

     if (mm_old /= mm) then   ! START OF NEW MONTH !!!!!

        call Code_timer(tim_before)

        !subroutines/data that must be updated every month

        call readdiss(newseason)

!rv3_5_14: aircraft from "newmonth"      if ( AIRNOX ) call aircraft_nox(newseason)

        if (MasterProc .and. DEBUG_UNI ) write(6,*) 'maaned og sesong', &
             numt,mm,mm_old,newseason,oldseason

        call Add_2timing(6,tim_after,tim_before,"readdiss, aircr_nox")

        call MetModel_LandUse(2)   ! e.g.  gets snow
        if ( MasterProc .and. DEBUG_UNI) write(6,*)"vnewmonth start"

        call newmonth

        call Add_2timing(7,tim_after,tim_before,"newmonth")

        if ( AIRNOX ) call lightning()

        call init_aqueous()


        call Add_2timing(8,tim_after,tim_before,"init_aqueous")

     end if    ! mm_old.ne.mm

     ! - we add a monthly call to BoundaryConditions. Can re-code later for
     !     possibly shorter call intervals

     call Code_timer(tim_before)

     if (mm_old /= mm) then   ! START OF NEW MONTH !!!!!
        if( DEBUG_UNI ) print *, "Into BCs" , me

        ! We set BCs using the specified iyr_trend
        !   which may or may not equal the meteorology year

        call BoundaryConditions(current_date%year,iyr_trend,mm)
        if( DEBUG_UNI ) print *, "Finished BCs" , me
        if(numt == 2) call Init_massbudget()
        if( DEBUG_UNI ) print *, "Finished Initmass" , me

     end if

     oldseason = newseason
     mm_old = mm

     call Add_2timing(9,tim_after,tim_before,"BoundaryConditions")

     if( DEBUG_UNI ) print *, "1st Infield" , me, " numu ", numt


     call Meteoread(numt)
        call Add_2timing(10,tim_after,tim_before,"Meteoread")

     call SetLandUse()    !Moved here from DryDep !LPJ
        call Add_2timing(11,tim_after,tim_before,"SetLanduse")

     call SetDailyBVOC(daynumber)  !dsBVOC

     if ( USE_FOREST_FIRES ) call Fire_Emis(daynumber)

        call Add_2timing(12,tim_after,tim_before,"Fires+BVOC")

     daynumber=day_of_year(current_date%year,current_date%month,&
          current_date%day)
     if ( MasterProc) write(6,*) 'TIME TEST ', 'current date ',current_date, &
             "day number ", daynumber


     call Code_timer(tim_before)

     call metvar(numt)

     call adv_var(numt)

     call Add_2timing(13,tim_after,tim_before,"metvar")

     call Code_timer(tim_before)

     call phyche(numt)

     call Add_2timing(14,tim_after,tim_before,"phyche")

     call WrtChem(numt)

     call trajectory_in


     call Add_2timing(37,tim_after,tim_before,"massbud,wrtchem,trajectory_in")


  end do !   end 3-hourly time-loop.

  call Code_timer(tim_after0)
  call Add_2timing(38,tim_after0,tim_before1,"total within loops")
  call Add_2timing(39,tim_after0,tim_before0,"total")

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !
  call wrtxn(current_date,.true.)

  call massbudget()

  if(MasterProc)then
     write(6,*) 'programme is finished'
    ! Gather timing info:
     if(NPROC-1 >  0)then
        CALL MPI_RECV( lastptim, 8*39, MPI_BYTE,  NPROC-1 &
             ,765, MPI_COMM_WORLD, STATUS, INFO)
     else
        lastptim(:) = mytimm(:)
     endif

     call Output_timing(IO_MYTIM,me,NPROC,nterm,GIMAX,GJMAX)

  else if(me == NPROC-1) then
     CALL MPI_SEND( mytimm, 8*39, MPI_BYTE,  0, 765, MPI_COMM_WORLD, INFO)
  endif
  !cccccccccccccccccccccccccccccccccc

  ! write 'modelrun.finished' file to flag the end of the FORECAST
  if (me==0 .and. FORECAST) then
    open(1,file='modelrun.finished')
    close(1)
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
  CALL MPI_FINALIZE(INFO)

end program myeul
