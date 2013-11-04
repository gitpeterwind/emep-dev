 module esx_Zchem
  use CheckStops, only : CheckStop
  use ChemSpecs, only : species, define_chemicals, nspec => NSPEC_TOT
  use ChemRates, only : NCHEMRATES
  use esx_Variables, only:  esx, Zmet
  use SmallUtils_ml, only: find_index, find_indices, NOT_FOUND
  use ZchemData,  only : Alloc1Dchem, xChem, rcemis

  implicit none
  private

  public  :: init_zchem
  private :: Config_Zchem

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0
   integer, private, parameter :: TXTLEN =   20

  !> Configuration of chemical boundary and initial conditions (BIC)
  !! =========================================================================
  !! Typically, we only need to specify BICs for long-lived species (advected). 

   integer, private, parameter :: MAXBICS = 20 !> Max. No. records with BICs

   type, private :: BIC_t
     character(len=10) :: type  = '-' ! ! conc or emis
     character(len=TXTLEN) :: name = "-"
     character(len=TXTLEN) :: units = '-' ! ! e.g. ppb, molec/cm2/s, but can be "-" for unitless
     real :: value= 0.0                !! concs in units given next:
     integer :: z1=0, z2=0             !! height level range for this record
   end type BIC_t

   ! Some more complex code for BIC moved to end (SKIP1)


  !! We store both types of initialisation for configZchem:

   type(BIC_t),    dimension(MAXBICS) :: BIC = BIC_t()


  contains
  !--------------------------------------------------------------------------!
  subroutine init_zchem(nz, ioconfig, debug_level) !!, errmsg)
    integer, intent(in) ::  nz   ! Number of layers
    integer, intent(in) :: ioconfig
    integer, intent(in) ::  debug_level  ! 0, 1 or 2
!    character(len=*), intent(inout) :: errmsg

    call Alloc1Dchem(nz, NCHEMRATES, debug_level)

    call Config_Zchem(ioconfig, debug_level ) !> Mainly boundary conditions

  end subroutine init_Zchem
  !--------------------------------------------------------------------------!
  !> The Config_Zchem routine is mainly to set initial and
  !! boundary conditions. We use text strings for compound names,

  subroutine Config_Zchem(ioconfig, debug_level)

    integer, intent(in) :: ioconfig
    integer, intent(in) ::  debug_level  ! 0, 1 or 2
    integer :: i, ilog=0, ispec, iz, iz1, iz2, nbic
    integer, dimension(MAXBICS) ::  bic_list = 0
    real :: ppb = 2.55e10    !! TMP fix later with real factor
    real :: unitscale        !! Converts input BICs to concs. 
    namelist /chem_config/BIC

    print *, "INTO ZCHEM==================================="
    rewind(ioconfig)
    read (ioconfig, nml=chem_config)
    if(  debug_level > 1 ) then
      open(newunit=ilog,file="LogConfig.Zchem")
      write(ilog,"(a)") "CONFIG chem ================================"
      write(ilog,nml=chem_config)
    end if

   !> Assign initial xChem values. Use simple approx for ppb just now.
   !!  Fix later !!
   !! ( Can use z2 = 999 to ensure top of domain )

    nbic = 0 

    BICLOOP : do i = 1,  MAXBICS

        print *, "BICLOOP ", i, trim( BIC(i)%name ), trim( BIC(i)%type )

        if ( BIC(i)%name == "-" ) exit
        ispec = find_index( BIC(i)%name, species(:)%name )
        if ( ispec < 1 ) then
          if(debug_level>0) write(*,*) "Zchem:BIC WARNING, BIC species not in mechanism:"&
                    //trim( BIC(i)%name)
          cycle BICLOOP
        end if

        if ( find_index( ispec, bic_list ) == NOT_FOUND) then ! new species
          nbic = nbic + 1
          if(debug_level>0)print *, "BICADDED ", i, trim( BIC(i)%name ), trim( BIC(i)%type ), nbic
          bic_list(nbic ) = ispec
        end if

        iz1 = BIC(i)%z1
        iz2 = min ( BIC(i)%z2, esx%nz )

        do iz = iz1, iz2

           if ( BIC(i)%units == "ppb" ) then
              unitscale = ppb !WILL MAKE z-dep later
           else if ( BIC(i)%units == "-" ) then
              unitscale = 1.0
           else if ( BIC(i)%units == "molec/cm2/s" ) then
              unitscale = 1.0/( 100.0 * esx%dz(iz) ) ! => /cm3/s
           else
              call CheckStop(.false., "ERROR: Wrong units in BIC:"//&
                         trim(BIC(i)%name)//":"//trim(BIC(i)%units) )
           end if
  
          if(debug_level>0) print "(2a,g12.3,2i4,es12.3)", "BIC unitscale:",&
             trim(BIC(i)%type), unitscale, iz1, iz2, BIC(i)%value
    
           if ( BIC(i)%type == 'conc' ) then
             xChem( ispec, iz ) = BIC(i)%value * unitscale
           else if ( BIC(i)%type == 'emis' ) then
             rcemis( ispec, iz ) = BIC(i)%value * unitscale
           else 

              call CheckStop(.false., "ERROR: Wrong type in BIC:"//&
                         trim(BIC(i)%name)//":"//trim(BIC(i)%type) )
           end if
!           if ( debug_level >1 ) write(ilog,"(a,3i4,es12.3,f8.2)") "BIC: "//trim(BIC(i)%type), &
!              i, ispec, iz, unitscale, BIC(i)%value
        end do

        if( debug_level >0 ) write(*, "(a,2i4,1x,a,2i3,g12.3)")  "BIC:SET ", &
           i, ispec, trim( BIC(i)%name ), iz1, iz2, BIC(i)%value

    end do BICLOOP

    if(  ilog /= 0 ) close(ilog)

    if( debug_level >0 ) then

      associate ( list => bic_list(1:nbic) )

        write(*,"(a,/,a,a3,2x,99a11)") &
          "ZChem_BIC Summary ---------------------",&
          "ZChem:BIC SET:", "iz", "C-"//species(list)%name, "E-"//species(list)%name
        do iz = esx%nz, 1, -1
            write(*,"(a,i3,99es11.3)") "Zchem:BIC SET:", iz, &
             xChem( list, iz), rcemis( list, iz)
        end do
      end associate
    end if


   ! Some more complex code for BIC moved to end (SKIP2)
     
  end subroutine Config_Zchem

  !--------------------------------------------------------------------------!
 end module esx_Zchem

!-----------------
!SKIP1
! Skip this fancier stuff for now. DiffSpecs now in esx_Variables
!  !> Configuration -  species affected by turbulent diffusion (ie by Kz)
!  !! =========================================================================
!  !! We have two  main options
!  !!   (a) "list" - specify a list of species names
!  !!   (b) "range" - specify a range of indices 
!  !! The latter is probably safer for large chemical schemes, and assumes that
!  !! the species in CM_ChemSpecs are ordered in some kind of lifetime scale. 
!  !! Typically we can skip all radicals, and also compounds like CO, CH4 etc.
!  !! (needs to be investigated)
!
!   integer, private, parameter :: MAXDIFF = 10 !> Max. No. species for diffusion
!   type, private :: KzSpec_t
!     character(len=TXTLEN)                   :: method = "list" !! list or range
!     character(len=TXTLEN), dimension(nspec) :: name = "-"
!     integer :: nesx = 0  !!  number of diffusing species for ESX runs
!     integer,  dimension(nspec) :: index = -1  !! index from species array
!   end type KzSpec_t
!-----------------
!SKIP2
!Skip this fancier stuff for now. DiffSpecs now in esx_Variables.
!   !> Calculate species for diffusion calcs.
!
!    associate ( list => Zchem%kzSpecs%name )
!    !list(:) = Zchem%kzSpecs%name(:) ! Copy list, and write back into KzSpecs%name

!    Zchem%KzSpecs%nesx = 0 !count( find_indices ( list(1:ispec), species(:)%name ) > 0 )

!    if ( Zchem%KzSpecs%method == "list" ) then

!      !sometimes the namelist will have species which are not found in the
!      ! current mechanism. We skip these

!      do i = 1, LenArray( list(:), "-" )   
!         ispec = find_index( list(i), species(:)%name  )
!         if ( ispec > 0 ) then
!           esx%nKzSpecs = esx%nKzSpecs + 1
!           Zchem%KzSpecs%index( esx%nKzSpecs ) = ispec
!DANGER            Zchem%KzSpecs%name( i ) = 
!         end if
!      end do
!
!    elseif ( Zchem%KzSpecs%method == "range" ) then

!      ispec1 = find_index( list(1), species(:)%name  )
!      ispec2 = find_index( list(2), species(:)%name  )
!      if( ispec2 < ispec1 ) call CheckStop(.true., "Zchem:KzSpecs:range NEG!")

!      do ispec = ispec1, ispec2
!           esx%nKzSpecs = esx%nKzSpecs + 1
!           Zchem%KzSpecs%index( esx%nKzSpecs ) = ispec
!           !DANGER Zchem%KzSpecs%name( esx%nKzSpecs ) = ispec
!      end do

!    else
!      call CheckStop(.true., "Zchem:KzSpecs:Incorrect method:"//&
!                     trim(Zchem%KzSpecs%method)) 
!    end if
!
!    do ispec = 1, esx%nKzSpecs
!       i = Zchem%KzSpecs%index( ispec )
!       write(*,*)"ESX KzSpecs ", ispec, i , trim(species(i)%name)
!    end do
!   end associate ! list
