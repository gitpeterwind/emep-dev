program tester
  use AllocInits,       only : AllocInit
  use Ammonium_ml,      only : Ammonium   ! Testing equilib system
  use ChemSpecs,        only : NSPEC_TOT, define_chemicals, species
  use CheckStop_ml,       only : CheckStop
  use esx_gleaf,        only : Set1Dgleaf
  use esx_GetData,      only : GetLocMetData    !! e.g. hourly surface met
  use esx_ChemRun,      only : ChemRun
  use esx_MassBudget,   only : mass, TOT, print_mass
  use esx_Variables,    only : Config_esx, esx, Zmet, Zveg, Loc
  use esx_Zchem,        only : init_zchem
  use esx_ZdiffSolver,  only : ZdiffSolver, nDiffSteps
  use esx_Zgrid,        only : init_zgrid
  use esx_Zmet,         only : set_esxZmet
  use esx_Zveg,         only : Config_Zveg, init_Zveg, Set1dPAR
  use Io_ml,            only : IO_LOG
  use Io_Progs,         only : PrintLog
  use Io_Routines,      only : writetdata
  use ZchemData,        only : Alloc1Dchem, xChem, rcemis
  use Zmet_ml,          only : Set1Dmet
  use Zmet_ml,          only : H2O ! for Plume2 tests
  use SmallUtils_ml,    only : find_index, trims, num2str
  use TimeDate_ml,      only : current_date, add2current_date

  implicit none

  character(len=500) :: plotmsg
  character(len=100) :: filename, sname
  character(len=50) :: txt
  real    :: units, pcm3_to_pm3, pm3_to_pcm3
  integer :: i,idspec, icspec, nprint, io_hour1, io_hour2, old_hour=-999
  logical :: first=.true.
  logical, save ::   my_first_call = .true.
  integer :: config_io
  integer :: nz, nDiffSpecs, nOutSpecs, nteval, nt

  !> Helper(or fake) arrays
  real, allocatable, dimension(:) :: cz  ! concs in 1-D column

  ! zeroth and 1st order source sink for Crank-Nichiolson, not used so far
  real, allocatable, dimension(:) :: SourceSink0, SourceSink1

  !> czprint Array stores results for printing/plotting
  real, allocatable, dimension(:,:,:) :: czprint  !(nz,nDiffSpecs,nteval)

  !> Check that compilation flags give higher precision
  call CheckStop(digits(1.0)<50, &
     "COMPILED WRONGLY: Need double precision, e.g. gfortran -fdefault-real-8")


  !=========================================
   call define_chemicals()       !! Sets species dims, names, indices etc
  !=========================================

  open (newunit=config_io, file="config_esx.nml")

   call Config_esx(config_io, writelog=.true.)   !> Sets esx%nz, etc.

   open (newunit=IO_LOG, file=trim(esx%odir)//"/Log.esxtester")

   call init_Zgrid(IO_LOG)                       !> esx%nz, esx%z, etc, ..

   call config_Zveg(config_io, writelog=.true.)  !> Sets Zveg%zbnd, etc.
 
   !> Reads config file, and allocates chemical arrays xChem, rcemis etc.

   call init_Zchem(esx%nz, config_io, debug_level=esx%debug_Zchem )

  close (config_io)

  !/ Some output files:
  txt=trim(esx%odir)// "/OutputHourly"
  open (newunit=io_hour1,file=trims( txt // "1m " // esx%exp_name // ".csv") )
  open (newunit=io_hour2,file=trims( txt // "top " // esx%exp_name // ".csv") )
!  open (newunit=io_hour2,file="OutputHourlytop" // trim(esx%exp_name) // ".csv")
  write (io_hour1,'(a,9999(:,",",a10))') "date", species(:)%name
  write (io_hour2,'(a,9999(:,",",a10))') "date", species(:)%name

  nz         = esx%nz
  nDiffSpecs = esx%nDiffSpecs
  nOutSpecs  = esx%nOutSpecs
  nteval = maxloc(esx%printTime,dim=1)    ! Number of print-out times
  write(txt,"(a,i4, ' zmax=',f8.2,' nteval=',i6)") "CONF Z", nz, esx%z(nz), nteval
  call PrintLog( txt )
  !write(*,"(a,i12,f15.3)") "TLOOPSTART", nDiffSteps

  call AllocInit( czprint,     0.0, nz,nOutSpecs,nteval, "czprint")
  call AllocInit( cz,          0.0, nz, "cz")
  call AllocInit( SourceSink0, 0.0, nz, "SS0")
  call AllocInit( SourceSink1, 0.0, nz, "SS1")

  call init_Zveg()    !> dLAI, cumLAI, => LogLAIz.txt output

  

!/ Solve 

  associate(t => esx%Time, tmax => esx%endTime, &
            dt=>esx%dt_phychem, tprint=>esx%printTime)

    TIMELOOP: do while (t+dt <= tmax)


      if( esx%debug_driver > 0 ) write(*,"(a,i12,f15.3)") "TLOOP", nDiffSteps, t

   !> Update external or large-scale met data for this time-step.

      call GetLocMetData(current_date)


   !> Update time-step met data for this time-step
   !! (Confusing with 2 calls I admit, but the first is meant to read in for
   !!  example EMEP model or experimental data, which ESX now has to convert
   !!  to 1-D arrays. The first is also similar to a routine in the EMEP
   !!  model. Will try to harmonise...)

      call Set_esxZmet()         !> basic met, temp (tzK), rhz, Kz, => Zmet
      call Set1Dmet( nz, Zmet )  !> Includes also derived arrays for chem,
                                 !!  e.g. O2, tinv

      txt =  "DMET P1,Pnz:" // num2str( Zmet(1)%Pa, '(es10.3)') // &
              num2str( Zmet(nz)%Pa, '(es10.3)')
      call PrintLog( txt )

      if( esx%exp_name == "Plume2" ) H2O = 2.55e18 !! TESTING

   !> Update veg PAR, gs etc if needed: (move inside time-loop soon)

      if( esx%uses_veg) then
          call Set1dPAR  ()
          call Set1Dgleaf() 
      end if
  !----------------- for printouts and logs -----------------------------------
   select case ( esx%units )
   case ( "ppb" ) 
     units = 1.0e9/Zmet(1)%M   !! ppb at surface
     pcm3_to_pm3 = 1.0e6   ! #/cm3 -> #/m3 inside KzSolver
   case ( "-" ) 
     units = 1.0
     pcm3_to_pm3 = 1.0
   case default
     stop "DEFINE esx%units"
   end select
   pm3_to_pcm3 = 1.0/pcm3_to_pm3
 
   if ( my_first_call ) then
     nprint = 1 
     do i = 1, nOutSpecs
       icspec = esx%OutSpecs(i)%int
       if( icspec < 1 ) then
           print *, "ICSPEC not found: ", icspec, esx%OutSpecs(i)%key
       else ! ok
           czprint(:, i , nprint) = xChem(icspec,:)*units
       end if
     end do
     nprint = nprint + 1

!/  Massbudget, initialise

     do i = 1, NSPEC_TOT
       mass(i)%init   = dot_product( xChem(i,:), esx%dz(1:esx%nz) )
     end do
     my_first_call = .false.
   end if ! my_first_call
   print *, "UNITS ppb ", Zmet(1)%M, units
!stop 'PPB'
  !----------------- for printouts and logs -----------------------------------

      !> CHEMISTRY ======================================:
      !> The chemical scheme includes any emission terms, typically
      !! the main cause of  fast production/loss in the ESX system
      !! (and also of strong vertical gradients)

      if( esx%uses_chem) then

              print "(a, es18.5)", "PRECHEM ",  xChem( 18, 5)
          call ChemRun( dt, esx%debug_Zchem )
              print "(a, es18.5)", "POSCHEM ",  xChem( 18, 5 )
          call Ammonium( esx%debug_Zchem>0 )
      end if

      !> DISPERSION =====================================:
      !! Here we include only those species specified in
      !! esx%DiffSpecs arrays

      if( esx%uses_diff) then

         DIFFSPECS: do idspec = 1, nDiffSpecs !> LOOP OVER DISPERSING SPECIES

           associate ( dspec => esx%DiffSpecs(idspec) ) 
              icspec = dspec%ind

            !! chem has molec/cm3.
            !! Needs to be in molec/m3 to match units of z, Kz, etc.

              cz(:) = xChem( icspec, : ) * pcm3_to_pm3

            !! TMP. We should have pollutant-specific gleaf values, but for now just do

              if ( esx%uses_veg ) then
                SourceSink1(:) = -Zveg(1:nz)%dLAI * Zveg(1:nz)%gleaf
              end if

         !! Emissions can be handled with Zdiff, but not if we use chem too
              if ( .not. esx%uses_chem ) then
                 SourceSink0(1:nz) = rcemis( icspec, 1:nz )
              end if
           
              call ZdiffSolver(nz, icspec, dt, &
                 Vd=dspec%Vd,  Ve=dspec%Ve, &
                 D=SourceSink1, E=SourceSink0, &
                 Fb=dspec%Fb, Ft=dspec%Ft, fixedBC=esx%fixedBC, &
                 nSubSteps=nt, concn=cz, debug_level=esx%debug_Zdiff )

              xChem(icspec,:) = cz(:) * pm3_to_pcm3

           end associate ! dspec

         end do DIFFSPECS
         nDiffSteps = nDiffSteps + nt  ! Keep track of diffusion steps after loop
      end if 


      t = t + dt     !! Update time-step
      call add2current_date(current_date, dt)

      if ( current_date%hour /= old_hour ) then
       !print *, "CDATE :", current_date
       write(txt,fmt="(i4,3i2.2)" ) current_date%year,current_date%month,&
                current_date%day,current_date%hour
       write (io_hour1,'(a,9999(:,",",es10.3))') trim(txt), xChem(:,1)*units
       write (io_hour2,'(a,9999(:,",",es10.3))') trim(txt), xChem(:,nz)*units
       old_hour=current_date%hour
      end if

      

      !! Store for print out if needed
      if ( t>0.0 .and. abs( t-tprint(nprint)) < 1.0e-3 ) then !! Store output

          do i = 1, nOutSpecs
            icspec = esx%OutSpecs(i)%int
            sname  = esx%OutSpecs(i)%key
            print *, trim("CZPRINT " // sname), xChem(icspec,1), units, esx%units
            print "(a20,2i3,f8.2,i4,es12.3,a,2es12.3)", "CZPRINT " // sname, &
              nprint,  i, t, icspec, xChem(icspec,1)*units, &
              "    MASS ", dot_product( xChem(icspec,:), esx%dz(1:esx%nz) ), &
                  mass(icspec)%b(TOT)
            czprint(:, i , nprint) = xChem(icspec,:)*units
          end do
          nprint = nprint + 1
      end if


      first = .false.


    end do TIMELOOP

!!print *, "TEST TRIMS ", trim(trims ("abc def")) // "XX"
!print *, "TEST TRIMS ", trim(trims (" abc def "))// "XX"
!stop 'TR'

    do i = 1, nOutSpecs !> Output to files

       icspec = esx%OutSpecs(i)%int
       sname  =trim(esx%OutSpecs(i)%key )

       call print_mass( icspec, nDiffSteps, nz, &
                         xChem(icspec,:), species(icspec)%name )

       txt =trims( esx%exp_name // "_" // sname  ) ! label for ascii and plot files
       filename=trims(esx%odir // "/Results_"// txt // ".txt" )

       call writetdata(filename, tprint, esx%z(1:nz), czprint( :,i, 1:nprint-1) )

       if ( esx%uses_plotting) then 
         write(plotmsg,"(9a)") trim( esx%plot_cmds ), " -i ", trim(filename),&
          " -c ", trim(sname), " -o ", trims(esx%odir // "/PlotRes_"// txt // ".png" )
         print *, "Plot cmds =", trim(plotmsg)
         call system(trim(plotmsg))
       end if
    end do
  end associate

end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
