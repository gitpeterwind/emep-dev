!*****************************************************************************!
!> MODULE GasParticleCoeffs_mod
!..............................................................................
! Specifies data for deposition modelling for specific gases and aerosols.
! as well as mapping arrays between real species and surrogates.
! All components in EMEP or ESX must be mapped to those in this module. (This
! is done using GenChem.)
! Data here are from:
! Wesely, 1989, Atmos. Environ., 23, No.6, pp. 1293-1304
! Massman, WJ, 1998, Atmos. Env., Vol. 32, 6, 1111-1127
! Equations from Seinfeld&Pandis, 1998
!
! (EMEP/ESX do not use the Wesely methodology, but makes some use of the H*,
!   f0 ideas.)
!..............................................................................

module GasParticleCoeffs_mod
  use AeroFunctions_mod,     only: DpgV2DpgN
  use CheckStop_mod,         only: CheckStop, checkValid, StopAll
  use ChemDims_mod
  use ChemSpecs_mod           ! speces and define_chemicals for self_test
  use Config_module,         only: MasterProc  ! emep's master
  use NumberConstants,       only: UNDEF_R
  use OwnDataTypes_mod,      only: TXTLEN_SHORT
  use PhysicalConstants_mod, only: PRANDTL, BOLTZMANN, FREEPATH, PI
  use SmallUtils_mod,        only: find_index
  implicit none
  private

! !/-- we define a type to map indices of species to be deposited
!   to the lesser number of species where Vg is calculated

 type, public :: dep_t
    character(len=TXTLEN_SHORT) :: &
         name        &! Species name
        ,surrogate    ! Surrogate species in calculated dep arrays
    real    :: setRate   ! if CDDEP_SET, give vg in m/s
  endtype dep_t

  include 'CM_DryDep.inc'  ! => CM_DDepMap = dep_t
  include 'CM_WetDep.inc'  ! => CM_WDepMap = dep_t


  public ::   GetDepMapping
  public ::   GasCoeffs
  public ::   ParticleCoeffs
  public ::   InitGasCoeffs
  public ::   InitParticleCoeffs
  public ::   WetCoeffs
  public ::   self_test

  real, private, parameter :: &
     DH2O    = 21.78e-6 &! m2/s at STP, Massman
    ,NU_AIR0 =  1.35e-5  ! Kin. viscosity of air, m2/s, 273K
  real, public, save       :: nu_air = UNDEF_R  ! Kin. viscosity of air, m2/s

  integer, public, parameter ::&
        NDRYDEP_GASES = 14     &! no. of gases in Wesely tables, DDdefs below
       ,NDRYDEP_AERO  = 11     &! no. of particles in DDdefs below
       ,NDRYDEP_DEF   = NDRYDEP_GASES + NDRYDEP_AERO ! gases + aerosol defs
     !mafor ,NDRYDEP_DEF   = 17      ! gases + aerosol defs ! MSK 26.01.2015 start

type, private :: DD_t
  character(len=16) :: name
  real :: Dx       ! diffusivity at STP, m2/s
  real :: DH2ODx   ! DH2O/Dx from Wesely Tables ! NOT USED YET
  real :: Hstar    ! M/atm    effective Henry coeff.
  real :: pe, K    ! from Wesely Tables, K in 1/(M s)
  real :: f0       ! Ozone reactivity indicator
  real :: Rm       ! mesophyll resistance (not yet used)
  ! Particle properties as in EMEP-ACP
  real :: umDpgV     ! Aerodynamic diameter for particles, volume, um
  real :: sigma    ! sigma of log-normal dist.  for particles
  real :: rho_p    ! particle density
  integer :: Gb     ! GerberClass ! 
end type DD_t

!Had:
! Mappings to DpgV types above, and Gerber types (see AeroFunctions).
! For Gerber (Gb), -1 indicates to use dry radius
!character(len=4), dimension(NSAREA_DEF) :: &
!  SLABELS=[character(len=4)::'SIAF','PMF','SSF','DUF','PM','SSC','DUC','ORIG']
!integer, dimension(NSAREA_DEF) ::&
!  Inddry = [ 1, 1, 1, 1, 2, 3, 4, 3], &
!  Gb     = [ 1, 1, 2,-1, 1, 2,-1,-1]
!end type aero_t
!type(aero_t), public, save :: AERO = aero_t()


type(DD_t), public, dimension(NDRYDEP_DEF), parameter :: DDdefs = [ &
! Dx values not from Massman are set using Wesely ratios for now. Will
! calculate later.
! DD_t( 'NO   ',18.02e-6 , 1.3, 2.0E-03, 9999, 1.0E-02,0,999.,-1,-1,-1,-1,& !No need for NO
! QUERY - what is K? not used so far, but check for H2O
! Rm from Zhang et al 2002. Note Zhang has 0 for MVK,100 for MACR
! DpgV, rhop are particle diameter and density, from EMEP
!               Dx (m2/s)  DH2O   H*       pe         K   f0  Rm umDpgV sig rhop Gb
!                          /Dx
!--------------------------------------------------------------------
! Gases:
  DD_t( 'SO2  ',10.89e-6 , 1.9, 1.0E+05,   -5, 1.0E+04,   0,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'O3   ',14.44e-6 , 1.6, 1.0E-02,   28, 6.0E+08,   1,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'NO2  ',13.61e-6 , 1.6, 1.0E-02, 9999, 2.0E+06, 0.1,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'H2O  ',DH2O     , 1.0, 1.0e+14, 9999,      -1,   0,999.,  -1,-1,-1,-1)&
 ,DD_t( 'HNO3 ',DH2O/1.9 , 1.9, 1.0E+14,    7, 1.0E-02,   0,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'H2O2 ',DH2O/1.4 , 1.4, 1.0E+05,   23, 7.0E+00,   1,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'ALD', DH2O/1.6 , 1.6, 1.5E+01,   -1, 1.0E+04,   0,100.,  -1,-1,-1,-1)& !Acetaldehyde
 ,DD_t( 'HCHO ',DH2O/1.3 , 1.3, 6.0E+03,   -3, 1.0E+04,   0,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'ROOH',DH2O/1.6 , 1.6, 2.4E+02, 9999, 2.0E+00, 0.1,100.,  -1,-1,-1,-1)& !Methyl hydroperoxide(was OP)
 ,DD_t( 'PAA  ',DH2O/2.0 ,   2, 5.4E+02, 9999, 6.0E+02, 0.1,  0.,  -1,-1,-1,-1)& !Peroxyacetic acid, 10
 ,DD_t( 'ORA',DH2O/1.6 , 1.6, 4.0E+06,   -8, 1.0E+04,   0,  0.,  -1,-1,-1,-1)& !Formic acid'
 ,DD_t( 'NH3  ',19.78e-6 ,   1, 1.0E+05, 9999, 1.0E+04,   0,  0.,  -1,-1,-1,-1)& 
 ,DD_t( 'PAN  ',DH2O/2.6 , 2.6, 3.6E+00, 9999, 3.0E+03, 0.1,  0.,  -1,-1,-1,-1)&
 ,DD_t( 'HNO2 ',DH2O/1.6 , 1.6, 1.0E+05,    6, 4.0E-04, 0.1,  0.,  -1,-1,-1,-1)&
! Particles:
!A2018 SHOULD CHECK and make consistent with ACP2012 Table 6 (or updated version)
!               Dx (m2/s)  DH2O   H*   pe    K   f0  Rm  umDpgV    sig  rhop Gb Inddry
 ,DD_t( 'PMfS'  ,UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.33,1.8, 1600,1)& ! as SAI_F 
 ,DD_t( 'PMfNO3',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.33,1.8, 1600,1)&! as SIA_F
 ,DD_t( 'PMfNH4',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.33,1.8, 1600,1)&! as SIA_F
 ,DD_t( 'SSf'   ,UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.33,1.8, 1600,2)& ! as SSF
 ,DD_t( 'DUf'   ,UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.33,1.8, 1600,-1)& ! NEW?? CHECK
!A2018 QUERY
 ,DD_t( 'PMc  ',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 3.00,2.0, 2200,1)& ! as PM QUERY 20
 ! SSc, DUc and POLLd have dummy values, CHECK!!
 ,DD_t( 'SSc  ',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 3.00,2.0, 2200,2)& ! or rho 2000?
 ,DD_t( 'DUc  ',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 3.00,2.0, 2200,-1)&! or 2.2,2600?
 ,DD_t( 'POLLd',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 3.00,2.0, 2200,-1)&
 ,DD_t( 'nuc  ',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.008,2.0, 1400,-1)&
 ,DD_t( 'ait  ',UNDEF_R  , -1,   -1,   -1,  -1,  -1, 0., 0.06,2.0, 1200,-1)&
]

! from EMEP Aqueous
!/ SUBCLFAC is A/FALLSPEED where A is 5.2 m3 kg-1 s-1,
!  and the fallspeed of the raindroplets is assumed to be 5 m/s.
real, parameter :: FALLSPEED = 5.0               ! m/s
real, parameter :: SUBCLFAC = 5.2 / FALLSPEED

!/ e is the scavenging efficiency (0.02 for fine particles, 0.4 for course)
real, parameter :: EFF25 = 0.02*SUBCLFAC, &
                   EFFCO = 0.4 *SUBCLFAC, &
                   EFFGI = 0.7 *SUBCLFAC

type, private :: WD_t
  character(len=16) :: name
  real :: W_sca
  real :: W_sub
end type WD_t

integer, parameter :: NWETDEP_DEF = 22
type(WD_t), public, dimension(NWETDEP_DEF),parameter :: WDdefs = [ &
  WD_t('SO2'  , 0.3,  0.15)  &! Berge+Jakobsen
 ,WD_t('SO4'  , 1.0,  EFF25) &! Berge+Jakobsen
 ,WD_t('NH3'  , 1.4,  0.5 )  &! subcloud = 1/3 of cloud for gases
 ,WD_t('HNO3' , 1.4,  0.5)   &!
 ,WD_t('H2O2' , 1.4,  0.5)   &!
 ,WD_t('HCHO' , 0.1,  0.03)  &!
 ,WD_t('ECfn' , 0.05, EFF25) &
 ,WD_t('SSf'  , 1.6,  EFF25) &
 ,WD_t('SSc'  , 1.6,  EFFCO) &
 ,WD_t('SSg'  , 1.6,  EFFGI) &
 ,WD_t('PMf'  , 1.0,  EFF25) &!!
 ,WD_t('PMc'  , 1.0,  EFFCO) &!!
 ,WD_t('POLLw', 1.0,  SUBCLFAC) &! pollen
!RB extras:
!perhaps too high for MeOOH? About an order of magnitude lower H* than HCHO:
 ,WD_t('ROOH' , 0.05, 0.015) &! assumed half of HCHO - 
 ,WD_t('0p2' , 0.2, 0.06) &!
 ,WD_t('0p3' , 0.3, 0.09) &!
 ,WD_t('0p4' , 0.4, 0.12) &!
 ,WD_t('0p5' , 0.5, 0.15) &!
 ,WD_t('0p6' , 0.6, 0.18) &!
 ,WD_t('0p7' , 0.7, 0.21) &!
 ,WD_t('0p8' , 0.8, 0.24) &! 
 ,WD_t('1p3' , 1.3, 0.39) &!
]


!integer, public, parameter :: CDDEP_SET = -99

! Mapping arrays:

  type, private :: depmap_t
    integer :: idef   ! index in definitions arrays
    character(len=TXTLEN_SHORT) :: name
    logical :: is_gas  ! true if gas
    integer, allocatable, dimension(:) :: advspecs  ! list of IXADV valuesadvect
  end type depmap_t

  type(depmap_t), allocatable, dimension(:), public ::&
    DDmapping, WDmapping, DepMapping

  integer, public, save :: nddep, nwdep ! will be number of DDspec, WDspec
  integer, public, save, dimension(NSPEC_TOT) :: itot2DDspec = 0
  integer, public, save :: idcmpHNO3, idcmpO3, idcmpNH3, idcmpSO2, idcmpNO2
  integer, public, save :: idcmpPMfS, idcmpPMfNO3, idcmpPMfNH4

  !/ Some terms needed for gases+aerosols. we just use one container for all
  type, private :: ddep_t
    character(len=TXTLEN_SHORT) :: name
    real :: &
        Dx       &! diffusivity at STP, m2/s
       ,DH2ODx   &! DH2O/Dx from Wesely Tables
       ,Hstar    &! M/atm    effective Henry coeff.
       ,f0       &! Ozone reactivity indicator
       ,DxDO3    &! Dx / DO3 ratio
       ,Re       &! Reynolds
       ,Schmidt  &! 
       ,Grashof  &! Grashofs
       ,Rb_cor   &! two-thirds power of the Schmidt to Prandtl numbers (for Rb)
       ,sigterm  &! for self-test and printouts
       ,sigma    &! sigma of log-normal dist.  for particles
       ,lnsig2   &! log(sigma)**2  used for settling velocity
       ,rho_p    &! particle density
       ,DpgV     &! volume-based geometric mean diameter (m)
       ,DpgN      ! number-based geometric mean diameter (m)
    integer :: Gb     ! GerberClass ! 
    logical :: is_gas 
  end type ddep_t  
  type(ddep_t), public, allocatable, dimension(:), save :: DDspec

  type, private :: wdep_t
    character(len=TXTLEN_SHORT) :: name
    real :: &
        W_sca    &! diffusivity at STP, m2/s
       ,W_sub     ! DH2O/Dx from Wesely Tables
    logical :: is_gas 
  end type wdep_t  
  type(wdep_t), public, allocatable, dimension(:), save :: WDspec


!/ Stored here for self-test and printouts
!  real, private, dimension(NDRYDEP_GASES+1: NDRYDEP_DEF), save :: sigterm, Dg
  character(len=100), private  :: fmt = "(a20,10(2x,a,es10.2))"

contains
 ! Mapping. Dry and Wet use very similar setups, but we need to grab
 ! info from different Tables. 

  subroutine GetDepMapping()

    integer :: iDep, i,n, iadv, itot, idef, nidef, nrow, irow
    character(len=TXTLEN_SHORT), allocatable, dimension(:) :: defnames
    character(len=TXTLEN_SHORT) :: defname, advname
    type(dep_t), allocatable, dimension(:) :: CM_DepMap   !Dry or wet
    !type(depmap_t), dimension(:), pointer :: ptrDepMapping
    character(len=4) :: dcase
    integer, dimension(900) :: ni                  ! tmp array
    integer, dimension(900,NSPEC_ADV) :: iDepUsed  ! tmp array
    character(len=*), parameter :: dtxt='GetDepMap:'


   ! we fill the CM_DepMap with the surrogate and then the
   ! list of  advected species associated with that.

    DRYWETLOOP: do idep = 1, 2  ! Dry then Wet

      if (idep==1) then
        dcase = 'dry:'
        allocate(defnames(size(DDdefs)), CM_DepMap(NDRYDEP_ADV))
        defnames = DDdefs(:)%name
        CM_DepMap  = CM_DDepMap
      else
        dcase = 'wet:'
        allocate(defnames(size(WDdefs)), CM_DepMap(NWETDEP_ADV))
        CM_DepMap  = CM_WDepMap
        defnames = WDdefs(:)%name
      end if

      n=0
      nrow = 0
      iDepUsed(:,:) = -999
      ni(:) = 0  ! number adv elements for each row
      nidef=0    ! number idef defined so far

     ! read CM_XXXdep.inc and build up tmp array with e.g.
     ! idef iadv ....
     ! YYY -999  not used, wil skip later
     ! O3   O3
     ! PAN  PAN MPAN
     ! ROOH CH3OOH C2H5OOH ..

      DEPLOOP: do  i = 1, size(CM_DepMap) ! from CM_DryDep - wanted species
        advname = CM_DepMap(i)%name
        defname = CM_DepMap(i)%surrogate ! must be in DDdefs, WDdefs
        itot    = find_index( advname, species(:)%name ) ! -ve if not unique
        iadv    = itot - NSPEC_SHL
        idef    = find_index( defname, defnames )

       ! Surrogate needs to exist in Defs list and species in species list
        if ( iadv <1 .or. idef <1 ) then
          print '(2(a,2i5))', dtxt//'NEG DepMap index'//dcase//&
            trim(advname), iadv, idef, trim(defname)
!print *, 'DEFNAMES ', defnames
          call StopAll(dtxt//'NEG DepMap index'//dcase//CM_DepMap(i)%name )
        end if

        ! And just a safety check, since setRate hasn't been re-implemented 
        ! yet in EMEP system

        call CheckStop( any(CM_DepMap(:)%setRate > 0.0 ), &
                        dtxt//dcase//"ATTEMPT at setRate ")

        if ( ni(idef) ==0 ) nrow = nrow + 1 ! new surrogate found

        ni(idef) = ni(idef) + 1
        iDepUsed(idef,ni(idef)) = iadv  ! adds species to row for nide

      end do DEPLOOP

      if( idep==1) allocate(DDmapping(nrow))
      if( idep==2) allocate(WDmapping(nrow))
      allocate(DepMapping(nrow))

      irow=0
      do idef = 1, len(defnames)
        if ( ni(idef) > 0 ) then
          irow = irow + 1 

          DepMapping(irow)%idef = idef
          DepMapping(irow)%name = defnames(idef)

          allocate( DepMapping(irow)%advspecs(ni(idef)) ) ! IXADV values, one or more
          if(idep==1) allocate( DDmapping(irow)%advspecs(ni(idef)) )
          if(idep==2) allocate( WDmapping(irow)%advspecs(ni(idef)) )
          DepMapping(irow)%advspecs = [ ( iDepUsed(idef,n), n=1,  ni(idef)) ]

          if(MasterProc) then
            if( irow==1) write(*,'(a)') dtxt//dcase// '  --- Mapping ----'
            write(*,'(a,2i4,1x,a8,a)',advance='no') dtxt//dcase, &
                irow, idef, defnames(idef),': '
          end if

          do n = 1,  ni(idef)
            itot = DepMapping(irow)%advspecs(n) + NSPEC_SHL
            if(MasterProc) write(*,'(a8)',advance='no') species(itot)%name
            !Chem2DDspec(itot) = irow  ! Mapping from 'real' to DDspec

            if( idep==1 ) then ! for dry dep we need reverse mapping
              itot2DDspec(itot) = irow  ! Mapping from 'real' to DDspec
            end if
          end do
          if(MasterProc) write(*,*)' '
        end if
    end do
       
  ! Now, copy DepMapping to wet or Dry
   if (idep==1) DDmapping = DepMapping
   if (idep==2) WDmapping = DepMapping

   deallocate(defnames)
   deallocate(CM_DepMap)
   deallocate(DepMapping)

 end do DRYWETLOOP

! Allocate other arrays used below
 nddep = size(DDmapping)
 nwdep = size(WDmapping)
 allocate(DDspec(nddep))
 allocate(WDspec(nwdep))
 DDspec(:)%name = DDmapping(:)%name ! assign here to help printouts
 WDspec(:)%name = WDmapping(:)%name
 !call GasCoeffs(298.15) ! Just sets DDspec names
 if(MasterProc) write(*,*) dtxt//'ALLOCS', size(DDspec), size(WDSPEC)

 ! Some species we use a lot:
 idcmpO3   = find_index('O3',DDspec(:)%name)
 idcmpHNO3 = find_index('HNO3',DDspec(:)%name)
 idcmpSO2  = find_index('SO2',DDspec(:)%name)
 idcmpNO2  = find_index('NO2',DDspec(:)%name)
 idcmpNH3  = find_index('NH3',DDspec(:)%name)
 idcmpPMfS   = find_index('PMfS',DDspec(:)%name)
 idcmpPMfNO3 = find_index('PMfNO3',DDspec(:)%name)
 idcmpPMfNH4 = find_index('PMfNH4',DDspec(:)%name)

  end subroutine GetDepMapping


!==========================================================
!>  GasCoeffs sets diffisivities etc for temperature
!!  Assumes near-surface for now. Will add rho later
!!  Massman suggest DH2O = 0.2178 (T/T0)**1.80 cm2/s where T0=273.15. 
!!  This T-dep fits changes in nu_air and kt from Garratt A3 very
!!  well also.

  subroutine InitGasCoeffs()
    integer :: iO3, icmp, idef
    real :: DxO3

     iO3 = find_index('O3',DDdefs(:)%name)
     DxO3 = DDdefs(iO3)%Dx

    ! copy across some values from DDdefs to shorter DDspec array

     do icmp = 1, size(DDmapping(:)%name )
       idef = DDmapping(icmp)%idef

       DDspec(icmp)%name  = DDdefs(idef)%name
       DDspec(icmp)%is_gas  = .false.

       if ( idef >  NDRYDEP_GASES ) CYCLE ! Just do gases here

       DDspec(icmp)%Dx      = DDdefs(idef)%Dx
       DDspec(icmp)%DxDO3   = DDdefs(idef)%Dx / DxO3
       DDspec(icmp)%Hstar   = DDdefs(idef)%Hstar
       DDspec(icmp)%f0      = DDdefs(idef)%f0
       DDspec(icmp)%is_gas  = .true.

       ! Sc and Rb_cor same for all T
       DDspec(icmp)%Schmidt = NU_AIR0 / DDdefs(icmp)%Dx
       DDspec(icmp)%Rb_cor  = (DDspec(icmp)%Schmidt/PRANDTL)**(2.0/3.0)

       if(MasterProc) write(*,fmt) "DD_Coeffs: "//DDspec(icmp)%name, &
       !print *, "DD_Coeffs: "//DDspec(icmp)%name, &
          "Dx=", DDspec(icmp)%Dx, "DxDO3=",DDspec(icmp)%DxDO3, &
          "Sc=", DDspec(icmp)%Schmidt, "Rb_cor=", DDspec(icmp)%Rb_cor

     end do

  end subroutine InitGasCoeffs

  subroutine GasCoeffs(T)
    real, intent(in)  :: T    ! /temp, K
  
    integer :: icmp, idef, iO3
    real :: Tcorr, DxO3

    Tcorr  = (T/273.15)**1.8

    do icmp = 1, size(DDmapping(:)%name )
      idef = DDmapping(icmp)%idef
      if ( idef >  NDRYDEP_GASES ) CYCLE ! Just do gases here
      DDspec(icmp)%Dx  = DDdefs(idef)%Dx * Tcorr
    enddo

    !real :: mu_air   ! Dynamic viscosity of air, kg/ms/s
    ! Jacobson, eqns 4.55,  for mu_a
    ! mu_air = (1.8325e-5)*(416.16/(T+120))*((T/296.16)**1.5)
    ! nu_air = mu_air/1.085  !QUERY CEH had 1.085- some standard rho?

    nu_air =  NU_AIR0 * Tcorr

  end subroutine GasCoeffs

  !----------------------------------------------------------------------------
  subroutine WetCoeffs()
    integer :: icmp, idef
    do icmp = 1, size(WDmapping(:)%name )
      idef = WDmapping(icmp)%idef
      WDspec(icmp)%name   = WDdefs(idef)%name
      WDspec(icmp)%W_sca  = WDdefs(idef)%W_sca
      WDspec(icmp)%W_sub  = WDdefs(idef)%W_sub
      if(MasterProc) write(*,fmt) "WD_Coeffs: "//WDspec(icmp)%name
    end do
  end subroutine WetCoeffs

  !----------------------------------------------------------------------------
  ! Particle diffusivity and Schmidt
  ! EMEP follows Binkowki+Shankar's approach, which integrates properties
  ! over a log-normal distribution. Eqns from B+S Appendix labelled
  ! below. 
  ! NB Notation is tricky B&S use Dp for diffusivity
  ! CAVEAT - so far this knows nothing about wet-radii, so we can calculate knut etc
  ! for known fixed dry radii

  subroutine InitParticleCoeffs(debug_flag)
    logical, intent(in), optional :: debug_flag
    logical, save :: debug = .false.
    real ::  knut, lnSig2
    integer :: idef, icmp
    if( present(debug_flag) ) debug = debug_flag

    do icmp = 1, size(DDmapping(:)%name )
      idef = DDmapping(icmp)%idef

      if ( idef <= NDRYDEP_GASES )  CYCLE ! Only particles here

      DDspec(icmp)%sigma = DDdefs(idef)%sigma         !just a copy
      DDspec(icmp)%rho_p = DDdefs(idef)%rho_p         !just a copy
      DDspec(icmp)%Gb    = DDdefs(idef)%Gb            !just a copy
      DDspec(icmp)%DpgV  = 1.0e-6 * DDdefs(idef)%umDpgV ! um -> m

     !... volume median diameter (Dp in EMEP notation. Will change!
     ! DpgN = Dpg in Seinfeld&Pandis
     ! -> geometric number median diameter ! S&P, 7.52:
     !TEST DpgV(i) = 1.0e-6 ! for comp to Seinfeld

      DDspec(icmp)%DpgN  = DpgV2DpgN(DDspec(icmp)%DpgV, DDspec(icmp)%sigma)

     !QUERY DDspec(icmp)%Dg = exp (log( DDdefs(idef)%Dp ) - 3*lnSig2 )
     ! done  above I think, we DpgN as result

      lnSig2 = log( DDspec(icmp)%sigma )**2
      DDspec(icmp)%lnsig2  = lnSig2

      !... mass median diameter -> geometric diameter 
            !CHECK Dp DDspec(icmp)%Dg = exp(log(DDspec(icmp)%Dp)-3* DDspec(icmp)%lnsig2)
      !DONE?  DDspec(icmp)%Dg = exp(log(DDspec(icmp)%DpgV)-3* DDspec(icmp)%lnsig2)
      if(MasterProc) write(*,fmt) "DD_Coeffs: "//DDspec(icmp)%name, &
              "Dp=", DDspec(icmp)%DpgV, "DpgN=",DDspec(icmp)%DpgN

      !QUERY knut = 2*FREEPATH/DDspec(icmp)%Dg   ! Knut's number
      knut = 2*FREEPATH/DDspec(icmp)%DpgN   ! Knut's number

      DDspec(icmp)%sigterm = exp(-2.5*lnSig2)+1.246*knut*exp(-4*lnSig2) !A29,dpk,k=3 

      if(MasterProc) write(*,"(a,2i4,2a8,10es10.2)") &
        "PMi,sig,Dp,Dg,Kn,sigterm: ",icmp, idef,  trim(DDspec(icmp)%name), &
         trim(DDdefs(idef)%name), DDspec(icmp)%sigma, DDspec(icmp)%DpgV,&
         DDspec(icmp)%DpgN, knut, DDspec(icmp)%sigterm
    end do !icmp

  end subroutine InitParticleCoeffs

  subroutine ParticleCoeffs(T,rho,debug_flag)
  real, intent(in)  :: T    ! temp, K
  real, intent(in)  :: rho  ! air density, kg/m3
  logical, intent(in), optional :: debug_flag
  real :: Dpg, knut, lnSig2
  integer :: idef, icmp
  logical, save :: first_call = .true., debug

  !------------------------------------------------------
  if( first_call ) then 
     if( present(debug_flag) ) debug = debug_flag
     call checkValid( nu_air, "nu_air needed" ) 
  end if ! first call 
  ! ------------------------ first call ----------------

  !... Diffusion coefficient for poly-disperse , A29, A30 from
  !    Binkowski & Shankar

   do icmp = 1, size(DDmapping)
      idef = DDmapping(icmp)%idef
      if ( idef <= NDRYDEP_GASES )  CYCLE ! Only particles here

      Dpg =BOLTZMANN*T/(3*PI * nu_air *rho * DDspec(icmp)%DpgN) ! A30

      DDspec(icmp)%Dx  = Dpg * DDspec(icmp)%sigterm   ! A29, Dpk 

      DDspec(icmp)%Schmidt = nu_air / DDspec(icmp)%Dx

      if(MasterProc .and. first_call ) then
        write(*,"(a,2f7.3,10es10.2)") "PM Coeffs,T,rho,nu,Dp,Dg", &
           T,rho,  nu_air, DDdefs(idef)%umDpgV, DDspec(icmp)%DpgN
        write(*,"(a,10es10.2)") "  Polydisperse Dx,Sc:", &
           DDspec(icmp)%Dx,DDspec(icmp)%Schmidt
        Dpg = Dpg*(1+2*1.246*FREEPATH/DDspec(icmp)%DpgN)  !just for check
          !FAILED write(*,"(a,10es10.2)") "  Monodisperse Dx,Sc:", DpgV,nu_air/DpgV 
      end if

   end do
   first_call = .false.

  end subroutine ParticleCoeffs
  !----------------------------------------------------------------------------
  subroutine self_test()
  integer :: iT,icmp,idef
  real :: T=293.15, rho=1.2

    call GetDepMapping()
    call InitGasCoeffs()
    call InitParticleCoeffs()
    call GasCoeffs(T)
    call ParticleCoeffs(T,rho,debug_flag=.true.) !QUERY VALUES v LOW
    call WetCoeffs()

    write(*,*) "Values at 298K, units m2/s for Dx."
    write(*,"(a,es9.2)") " nu_air=", nu_air
    write(*,"(a,2i4)") "No gas and gas+particles ", NDRYDEP_GASES, NDRYDEP_DEF

    write(*,*) '-------------'
    do icmp = 1, size(DDmapping)
      idef = DDmapping(icmp)%idef
      if ( DDspec(icmp)%is_gas ) then
        write(*,fmt) "GasCoeffs "//DDdefs(idef)%name, "Dx=", DDspec(icmp)%Dx, &
           "Sc=", DDspec(icmp)%Schmidt, "Rb_cor=", DDspec(icmp)%Rb_cor, &
           "Zhang Rm=", DDdefs(idef)%Rm,  & ! NEEDS CHECK for Rm
           " W&Th. Rm? ", 1.0/(DDdefs(idef)%Hstar/3000.0 + 100*DDdefs(idef)%f0 )
      else ! particle
        write(*,fmt) "PmCoeffs "//DDdefs(idef)%name, "Dx=",DDspec(icmp)%Dx, &
         "DpgV=", DDspec(icmp)%DpgV, "Sig=", DDspec(icmp)%sigma, &
         "DpgN=", DDspec(icmp)%DpgN, "Sc=",DDspec(icmp)%Schmidt
      end if
    end do !icmp

    icmp=find_index('O3',DDdefs(:)%name )
    do iT = 0, 30, 5    
      T=273.15 + iT
      call GasCoeffs(T)
!TMP      write(*,fmt) DDdefs(n)%name//"=f(T) Coeffs ", "T=", T, "Dx ", Dx(n)
    end do !iT
  end subroutine self_test
end module GasParticleCoeffs_mod
!TSTEMXprogram testdd
!TSTEMX  use GasParticleCoeffs_mod, only : self_test
!TSTEMX  use ChemSpecs_mod           ! speces and define_chemicals for self_test
!TSTEMX  implicit none
!TSTEMX  call define_chemicals()
!TSTEMX  call self_test()
!TSTEMXend program testdd
