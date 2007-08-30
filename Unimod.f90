program myeul
  !
  !     this is the main program for the off-line regional scale multilayer
  !     eulerian model at emep/msc-w. the main program contains the outer
  !     time-loop which runs through all time-levels a new meteorological
  !     data-set is read into the model from file. the inner time-loop
  !     runs through the physical time-step.
  !
  !
  !     the structur of the model is the following:
  !
  !------------------------------------------------------------------------

  !
  !
  !              (1)  a short description of each subroutine. a more
  !                   detailed description is found in each subroutine.
  !
  !              (2)  declaration and explanation of all model parameters are
  !                   given in the include-files: eulpar.inc, eulcon.inc,
  !                   eulbud.inc, eulmc.inc
  !
  !              (3)  the first meteorological data- and parameterset is read
  !                   from file. other initial parameters and model fields are
  !                   also defined under this paragraph.
  !
  !
  !     ----------------------
  !     begin outer time-loop:
  !     ----------------------
  !
  !
  !              (4)  new meteorological fields are derived and all the
  !                   necessary interpolation in space and time is performed
  !                   here.
  !
  !           ----------------------
  !           begin inner time-loop:
  !           ----------------------
  !
  !
  !              (5)  the chemical and physical calculations are performed 
  !                   in this paragraph.
  !
  !           --------------------
  !           end inner time loop.
  use My_Emis_ml,       only : NBVOC, AIRNOX
  use My_Outputs_ml,    only : set_output_defs
  use My_Timing_ml,     only : lastptim,mytimm,Output_timing, &
       Init_timing, Add_2timing, Code_timer, &
       tim_before,tim_before0,tim_before1, &
       tim_after,tim_after0
  use My_WetDep_ml,     only : Init_WetDep
  use MyChem_ml,        only : Init_mychem   

  use Advection_ml,     only : vgrid,adv_var, assign_nmax,assign_dtadvec
  use Aqueous_ml,       only : init_aqueous   !  Initialises & tabulates
  use AirEmis_ml,       only : aircraft_nox, lightning
  use Biogenics_ml,     only : Init_BVOC
  use BoundaryConditions_ml, only : BoundaryConditions
  use CheckStop_ml,     only : CheckStop
  use Chemfields_ml  , only : xn_adv
  use DefPhotolysis_ml, only : readdiss
  use Derived_ml,    only :  Init_Derived &
       ,IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY
  use Emissions_ml,     only : Emissions ,newmonth      !  subroutines
  use GenChemicals_ml,  only : define_chemicals
  use GenSpec_adv_ml  , only : NSPEC_ADV
  use GridValues_ml,    only : MIN_ADVGRIDS,GRIDWIDTH_M,Poles
  use Io_ml  ,          only : IO_MYTIM,IO_RES,IO_LOG,IO_TMP
  use Io_Progs_ml  ,    only : read_line
  use Landuse_ml,       only : ReadLandUse
  use MassBudget_ml,    only : Init_massbudget,massbudget
  use Met_ml  ,         only : metvar,MetModel_LandUse,&
       tiphys,Meteoread,MeteoGridRead,&
       startdate
  use ModelConstants_ml,only : KMAX_MID   &   ! No. vertical layers
       ,NPROC      &   ! No. processors
       ,METSTEP    &   ! Hours between met input
       ,runlabel1  &   ! explanatory text
       ,runlabel2  &   ! explanatory text
       ,nprint,nass,nterm,iyr_trend, PT
  use NetCDF_ml,        only : InitnetCDF,Init_new_netCDF
  use out_restri_ml,    only : set_outrestri
  use OutputChem_ml,    only : WrtChem
  use Par_ml,           only : me,GIMAX,GJMAX ,MSG_MAIN1,MSG_MAIN2&
       ,Topology, parinit & 
       ,li0,li1,lj0,lj1,tgi0,tgj0  &   ! FOR TESTING
       ,limax, ljmax, MAXLIMAX, MAXLJMAX, gi0, gj0
  use PhyChem_ml,       only : phyche  ! Calls phys/chem routines each dt_advec
  use TimeDate_ml,      only : date, current_date, day_of_year,daynumber
  use Trajectory_ml,    only : trajectory_init,trajectory_in ! For extraction of
  ! data along 3-D trajectories

  use Sites_ml,         only : sitesdef  ! to get output sites
  use Tabulations_ml,   only : tabulate
  use Nest_ml,         only : wrtxn


  !           --------------------
  !
  !
  !              (6)  output of chemical and meteorological fields, averaging,
  !                   accumumlation and massbudgets.
  !
  !
  !     --------------------
  !     end outer time-loop.
  !     --------------------
  !
  !
  !---------------------------------------------------------------------------

  !
  !
  ! ++(1)   the subroutines are:
  !
  !
  !      Before the time loops the model topography is set, initial fields 
  !      are read etc.
  !
  !     parinit           - Find the x-, y-, and z-addresses of the domain 
  !                         assigned to the processor.
  !     set_outrestri     - (in Out_restri_ml), Defines and checks the
  !                         restricted output dommain
  !     define_chemicals  - sets up species details
  !
  !     set_output_defs   - Initialises outputs
  !
  !     assign_nmax       - Assign nr. of time-steps for the inner time-loop
  !
  !     trajectory_init      - Should be replaced by a more generad 3D-trajectory
  !                         (presently not used)
  !     Meteoread         - Reads in first set of meteorological data
  !
  !     Emissions         - reads the emission totals in each grid-square
  !
  !     metvar            - derive new meteorological quantities from the
  !                         met. data read from the field-files.
  !     tabulate          - tabulates rct to rctit, sets up tab_esat, etc.  
  !
  !     sitesdef          - see if any output for specific sites is wanted 
  !                        (read input files "sites.dat" and "sondes.dat")
  !     vgrid             - Inclusion of the variable grid spacing
  !
  !     metvar            - Derived meteorological fields calculated from basic fields
  !
  !     adv_var           - Exchange u,v info for shadow area on nodes etc??
  !
  !
  !     3-hourly time loop
  !     --------------------
  !
  !         if ( new month) then
  !
  !             readdiss     - Read in new sets of j-values every month
  !
  !             aircraft_nox - (optional) seasonal aircraft emissions
  !
  !             MetModel_LandUse  - Reads monthly snow cover and yearly roughness data
  !                                 as used in HIRLAM/xx met model.
  !
  !             newmonth     - Reads natural so2 emissions (Should the name be
  !                             changed to ie monthly_DMS ??? )
  !
  !             lightning    - (optional) reads lightning emissions.
  !
  !             init_aqueous - set up of aqueous rates (moved before loop?? )
  !
  !             tstfld       - now only for debugging (see below)
  !
  !         end if ( new month )
  !
  !         if ( new month) then
  !
  !             BoundaryConditions  - Update boundary concentrations (and
  !                                   possible background fields)
  !             if ( first time in loop )
  !                init_massbudget  - Preset to zero and calc. total
  !                                   initial mass
  !             end if ( first time in loop )
  !
  !         end (new month )
  !
  !         remaining routines are called for every nterm
  !
  !         Meteoread - Reads in a new set of meteorological data
  !
  !         metvar    - Derived fields calculated from basic fields
  !                     (from infield)
  !
  !         adv_var   - Exchange u,v info for shadow area on nodes etc??
  !
  !         phyche:   - performes the physical and chemical calculations,
  !                     an explanation is found in each routine.
  !
  !                     subroutines: - for setting daynumber, emissions,
  !                                    zenith angle
  !
  !                                  - Advection rutines ( advecdiff )
  !
  !                                  - runchem calls setup subs and runs
  !                                                            chemistry
  !                                  - drydep
  !
  !                                  - ++ ( see phyche )
  !
  !
  !         Wrtchem    - Writes to output files (range:  3 hourly - annual)
  !
  !         trajectory_in - See above
  !
  !     end 3-hourly time loop
  !     -------------------------
  !
  !     massbud:       - calculates and print out budget quantities
  !
  !     output_timing  - writes the CPU spent in various parts of the code 
  !
  !
  !--------------------------------------------------------------------
  !
  !
  !++(2) explanation of model parameters/arrays in alphabetic order. help
  !      variables of minor importance are not explained here.
  !
  !
  !     an             - number og grid-distances from equator to pole
  !     xp             - x coordinate of the pole
  !     yp             - y coordinate of the pole
  !     cp             - heat capacity of air at constant pressure
  !     dt_advec       - length of advection (phyche) time-step
  !     efac           - conversion factor for emission totals
  !                      (see subroutine emis)
  !     epsil          - epsilon (a very small value)
  !     GRIDWIDTH_M    - grid-distance
  !     i              - x directed index
  !     j              - y directed index
  !     k              - vertical coordinate index
  !     NSPEC_ADV      - number of chemical components transport in the model
  !     NSPEC_SHL      - number of chemical components not-transport in the model
  !     NSPEC_BGN      - number of specified chemical components
  !     r              - gas constant
  !     t1             - start real-time of model calculation
  !     t2             - end real-time of calculation
  !     ndays          - number of days since 1 january (max 365 or 366)
  !     thour          - utc-time in hours every time-step
  !     xm(i,j)        - map factor
  !     xm2(i,j)       - xm**2
  !     xmd(i,j)       - 1./xm2
  !     xn_adv(NSPEC_ADV,i,j,k) - chemical component (1 = so2, 2 = so4)
  !     .
  !     .
  !     .
  !     .
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

  integer d


  logical, parameter :: DEBUG_UNI = .false. 
  integer n, numt, nadd, oldseason,newseason
  integer iupw, i, j, ii, k, iotyp
  integer :: mm, mm_old   ! month and old-month  (was nn and nold)
  integer :: nproc_mpi,cyclicgrid
  character (len=130) :: fileName, errmsg,txt

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
     write(unit=errmsg,fmt=*)"Wrong processor number!", &
          "program was compiled with NPROC = ",NPROC, &
          " but linked with ", nproc_mpi 
     call CheckStop( errmsg )
  end if
  call CheckStop( digits(1.0) < 50, &
       "COMPILED WRONGLY: Need double precision, e.g. f90 -r8")

  if (me  ==  0) then
     open(IO_RES,file='eulmod.res')
     open(IO_LOG,file='RunLog.out')
     open(IO_TMP,file='INPUT.PARA')
  endif

  call read_line(IO_TMP,txt,status(1))
  read(txt,*) nterm
  call read_line(IO_TMP,txt,status(1))
  read(txt,*) nass
  call read_line(IO_TMP,txt,status(1))
  read(txt,*) iyr_trend
 
  call read_line(IO_TMP,runlabel1,status(1))! explanation text short
  call read_line(IO_TMP,runlabel2,status(1))! explanation text long
  do i=1,3
     call read_line(IO_TMP,txt,status(1))
     read(txt,*)startdate(i)! meteo year,month,day to start the run
  enddo
  startdate(4)=0!meteo hour to start the run

  if( me == 0 ) then
     close(IO_TMP)
     write(unit=IO_LOG,fmt=*)trim(runlabel1)
     write(unit=IO_LOG,fmt=*)trim(runlabel2)
     write(unit=IO_LOG,fmt=*)startdate(1)
     write(unit=IO_LOG,fmt=*)startdate(2)
     write(unit=IO_LOG,fmt=*)startdate(3)
     write(unit=IO_LOG,fmt=*)"iyr_trend= ", iyr_trend
  endif

  print *, "read standard input"
  if( me == 0 ) print *, "RUNLABEL INPUT ", trim(runlabel1),' ',trim(runlabel2)
  if( me == 0 ) print *, " Trend Year is ", iyr_trend


  !*** Timing ********

  call Init_timing()
  call Code_Timer(tim_before0)
  tim_before = tim_before0

  call parinit(MIN_ADVGRIDS)     !define subdomains sizes and position
  call MeteoGridRead(cyclicgrid) !define grid projection and parameters
  call Topology(cyclicgrid,Poles)      !define GlobalBoundaries 
  !and subdomains neighbors
  call assign_dtadvec(GRIDWIDTH_M)! set dt_advec

  call set_outrestri  ! Steffen's routine for restricted area output


  !     Decide the frequency of print-out
  !
  nadd = 0
  nprint = nterm
  if (nterm >  nprint) nadd = 1

  if (me ==  0) write(6,*)'nterm, nprint',nterm, nprint

  !-------------------------------------------------------------------
  !
  !++(3)  parameters and initial fields.
  !

  call Add_2timing(1,tim_after,tim_before,"Before define_Chemicals")

  call define_chemicals()    ! sets up species details

  call Init_Derived()        ! Derived field defs., rv1_9_16

  call set_output_defs()     ! Initialises outputs

  call assign_nmax(METSTEP)    !  nmax=????

  call trajectory_init

  call Add_2timing(2,tim_after,tim_before,"After define_Chems, readpar")


  call MeteoRead(1)

  call Add_2timing(3,tim_after,tim_before,"After infield")

  if ( me == 0 ) write(6,*)"Calling emissions with year" ,current_date%year
  call Emissions(current_date%year)  !!!!! IS this the right/best year????
  if ( me == 0 ) write(6,*)"emissions fifniseh" 

  ! daynumber needed  for BCs, so call here to get initial

  daynumber=day_of_year(current_date%year,current_date%month,current_date%day)


  call MetModel_LandUse(1)   !

  call ReadLandUse()

  if ( NBVOC > 0  ) call Init_BVOC()

  call tabulate()    ! =>  sets up tab_esat, etc.

  call Init_mychem()   ! tabulates rct to rctit

  call Init_WetDep()   ! sets up scavenging ratios



  call sitesdef()      !--- see if any output for specific sites is wanted
  !   (read input files "sites.dat" and "sondes.dat" )


  call vgrid    !  initialisation of constants used in vertical advection
  if ( me == 0 ) write(6,*)"vgrid fifniseh" 

  if ( me == 0 ) then
     fileName=trim(runlabel1)//'_inst.nc'
     iotyp=IOU_INST
     call Init_new_netCDF(fileName,iotyp) 
     !netCDF hourly is initiated in Output_hourly
     fileName=trim(runlabel1)//'_hour.nc'
     iotyp=IOU_HOUR
     call Init_new_netCDF(fileName,iotyp) 
     fileName=trim(runlabel1)//'_day.nc'
     iotyp=IOU_DAY
     call Init_new_netCDF(fileName,iotyp) 
     fileName=trim(runlabel1)//'_month.nc'
     iotyp=IOU_MON
     call Init_new_netCDF(fileName,iotyp) 
     fileName=trim(runlabel1)//'_year.nc'
     iotyp=IOU_YEAR
     call Init_new_netCDF(fileName,iotyp) 

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

        if ( AIRNOX ) call aircraft_nox(newseason)

        if (me == 0) write(6,*) 'maaned og sesong', &
             numt,mm,mm_old,newseason,oldseason

        call Add_2timing(6,tim_after,tim_before,"readdiss, aircr_nox")

        call MetModel_LandUse(2)   ! e.g.  gets snow
        if ( me == 0 ) write(6,*)"vnewmonth start" 

        call newmonth

        call Add_2timing(7,tim_after,tim_before,"newmonth")

        if ( AIRNOX ) call lightning()

        call init_aqueous()

        if(numt == 2) call tstfld  !test output


        call Add_2timing(9,tim_after,tim_before,"init_aqueous")

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

     call Add_2timing(8,tim_after,tim_before,"BoundaryConditions")

     if( DEBUG_UNI ) print *, "1st Infield" , me, " numu ", numt


     call Meteoread(numt)

     call Add_2timing(10,tim_after,tim_before,"infield")

     daynumber=day_of_year(current_date%year,current_date%month,&
          current_date%day)
     if ( me == 0) then
        write(6,*) 'TIME TEST ', 'current date ',current_date, &
             "day number ", daynumber
     endif


     !++(5)
     call Code_timer(tim_before)

     call metvar(numt)

     call adv_var(numt)

     call Add_2timing(11,tim_after,tim_before,"metvar")

     call Code_timer(tim_before)

     call phyche(numt)

     call Add_2timing(18,tim_after,tim_before,"phyche")

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

  if(me == 0)then
     write(6,*) 'programmet er ferdig'
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

  CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
  CALL MPI_FINALIZE(INFO)

  stop
end program myeul
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!
!
!     +------------------------------------------------------------+
!     +                                                            +
!     +                  end main program                          +
!     +                                                            +
!     +------------------------------------------------------------+
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine tstfld
  !     - assignes an initial test-distribution to a
  !       chemical component. Not currently used except as print-out
  !
  integer i, j, k, n, info
  real rwork
  logical, parameter :: DEBUG_TEST = .false. 

  if ( DEBUG_TEST ) then
     write(6,*) "TSTFIELD OUTPUTS: me ", me
     write(6,*) "MAXS:", me, MAXLIMAX, MAXLJMAX
     write(6,*) "GMAXS:", me, GIMAX, GJMAX
     write(6,*) "limaxs:", me, limax, ljmax
     write(6,*) "li0:", me, li0, lj0
     write(6,*) "li0:", me, li1, lj1
     write(6,*) "gi0:", me, gi0, gj0
     write(6,*) "tgi0:", me, tgi0(me), tgj0(me)
     write(6,*) "DONE WITH TESTFIELD FOR ME ", me
  end if

  ! Not used now. Kept for future tests.
  !    do k=1,KMAX_MID
  !      do j=lj0,lj1
  !        do i=li0,li1
  !            xn_adv(IXADV_TEST,i,j,k) = ??????
  !        enddo
  !      enddo
  !    enddo

end subroutine tstfld
