program tester
  use AllocInits,       only : AllocInit
  use Ammonium_ml,      only : Ammonium   ! Testing equilib system
  use ChemSpecs,        only : NSPEC_TOT, define_chemicals, species
  use CheckStops,       only : CheckStop
  use esx_gleaf,        only : config_DO3SE, DO3SE_conf => conf, Set1Dgleaf
  use esx_GetData,      only : GetExternData    !! e.g. hourly surface met
  use esx_ChemRun,      only : ChemRun
  use esx_MassBudget,   only : mass, TOT, print_mass
  use esx_Variables,    only : Config_esx, esx, Zmet, Zveg, Loc
  use esx_Zchem,        only : init_zchem
  use esx_ZdiffSolver,  only : ZdiffSolver, nDiffSteps
  use esx_Zgrid,        only : init_zgrid
  use esx_Zmet,         only : set_esxZmet
  use esx_Zveg,         only : Config_Zveg, init_Zveg, Set1dPAR
  use Io_ml,            only : IO_LOG
  use Io_Routines,      only : writetdata
  use ZchemData,        only : Alloc1Dchem, xChem, rcemis
  use Zmet_ml,          only : Set1Dmet
  use SmallUtils,       only : find_index
  use TimeDate,         only : current_date, add2current_date

  implicit none

  character(len=100) :: msg = "ok", filename, sname
  real    :: units
  integer :: i,idspec, icspec, nprint
  logical :: first=.true.
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

   call init_Zgrid()                             !> esx%nz, esx%z, etc, ..

   call config_Zveg(config_io, writelog=.true.)  !> Sets Zveg%zbnd, etc.
 
   !> Reads config file, and allocates chemical arrays xChem, rcemis etc.

   !call Alloc1Dchem(esx%nz, esx%debug_Zchem )

   call init_Zchem(esx%nz, config_io, debug_level=esx%debug_Zchem )

  close (config_io)

  !/ Some output files:
  open (newunit=IO_LOG, file="Log.esxtester")

  nz         = esx%nz
  nDiffSpecs = esx%nDiffSpecs
  nOutSpecs  = esx%nOutSpecs
  nteval = maxloc(esx%printTime,dim=1)    ! Number of print-out times
  print "(a,i4, ' zmax=',f8.2,' nteval=',i6)", "CONF Z", nz, esx%z(nz), nteval
  write(*,"(a,i5,f15.3)") "TLOOPSTART", nDiffSteps

  call AllocInit( czprint,     0.0, nz,nOutSpecs,nteval, "czprint")
  call AllocInit( cz,          0.0, nz, "cz")
  call AllocInit( SourceSink0, 0.0, nz, "SS0")
  call AllocInit( SourceSink1, 0.0, nz, "SS1")

  call init_Zveg()    !> dLAI, cumLAI, => LogLAIz.txt output

  ! testing
  
   nprint = 1 
   select case ( esx%units )
   case ( "ppb" ) 
     units = 1.0e9/Zmet(1)%M   !! ppb at surface
   case ( "-" ) 
     units = 1.0
   case default
     stop "DEFINE esx%units"
   end select

   do i = 1, nOutSpecs
     icspec = esx%OutSpecs(i)%int
     if( icspec < 1 ) print *, "ICSPEC NEG ", icspec, esx%OutSpecs(i)%key
     czprint(:, i , nprint) = xChem(icspec,:)*units
   end do
   nprint = nprint + 1

!/  Massbudget, initialise

   do i = 1, NSPEC_TOT
     mass(i)%init   = dot_product( xChem(i,:), esx%dz(1:esx%nz) )
   end do

!/ Solve 

  associate(t => esx%Time, tmax => esx%endTime, &
            dt=>esx%dt_phychem, tprint=>esx%printTime)

    TIMELOOP: do while (t+dt <= tmax)


      if( esx%debug_driver > 0 ) write(*,"(a,i5,f15.3)") "TLOOP", nDiffSteps, t

   !> Update external met data for this time-step. Assumed hourly for now

      if( esx%uses_ExternData) call GetExternData(current_date)


   !> Update time-step met data for this time-step
   !! (Confusing with 2 calls I admit, but the first is meant to be ESX 
   !!  specific, and the 2nd EMEP model consistent. Will try to harmonise...)

      call Set_esxZmet()         !> basic met, temp (tzK), rhz, Kz, => Zmet
      call Set1Dmet( nz, Zmet )  !> Includes also derived arrays for chem,
                                 !!  e.g. O2, tinv

   !> Update veg PAR, gs etc if needed: (move inside time-loop soon)

      if( esx%uses_veg) then
          call Set1dPAR  ()
          call Set1Dgleaf() 
      end if

      !> CHEMISTRY ======================================:
      !> The chemical scheme includes any emission terms, typically
      !! the main cause of  fast production/loss in the ESX system
      !! (and also of strong vertical gradients)

      if( esx%uses_chem) then

          call ChemRun( dt, esx%debug_Zchem )
          call Ammonium( esx%debug_Zchem>0 )
      end if

      !> DISPERSION =====================================:
      !! Here we include only those species specified in
      !! esx%DiffSpecs arrays

      if( esx%uses_diff) then

         DIFFSPECS: do idspec = 1, nDiffSpecs !> LOOP OVER DISPERSING SPECIES

           associate ( dspec => esx%DiffSpecs(idspec) ) 
              icspec = dspec%ind

              cz(:) = xChem( icspec, : )

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

              xChem(icspec,:) = cz(:)

           end associate ! dspec

         end do DIFFSPECS
         nDiffSteps = nDiffSteps + nt  ! Keep track of diffusion steps after loop
      end if 


      t = t + dt     !! Update time-step
      call add2current_date(current_date, dt)
      print *, "CDATE :", current_date
      

      !! Store for print out if needed
      if ( t>0.0 .and. abs( t-tprint(nprint)) < 1.0e-3 ) then !! Store output

          do i = 1, nOutSpecs
            icspec = esx%OutSpecs(i)%int
            sname  = esx%OutSpecs(i)%key
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


    do i = 1, nOutSpecs !> Output to files

       icspec = esx%OutSpecs(i)%int
       sname  =trim(esx%OutSpecs(i)%key )

       call print_mass( icspec, nDiffSteps, nz, &
                         xChem(icspec,:), species(icspec)%name )

       filename="Results_"//trim(esx%exp_name) // "_" // trim(sname ) // ".txt"
       call writetdata(filename,  tprint, esx%z(1:nz),  czprint( :,i, 1:nprint-1) )

       if ( esx%uses_plotting) then 
         write(msg,"(9a)") trim( esx%plot_cmds ), " -i ", trim(filename), " -c ", trim(sname)
         print *, "Plot cmds =", trim(msg)
         !call system("cat "//trim(filename) )
         call system(trim(msg))
       end if
    end do
  end associate

end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
