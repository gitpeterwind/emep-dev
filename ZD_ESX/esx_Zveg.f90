!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> MODULE esx_Zveg
!>   Definition of leaf area variables: nhVeg, dLAI(z), cumLAI(z)
!!   gsto, gns, PAR.
!!   LAI - conversion from JP matlab routines, May 2013. (minor changes)
!@author
!> David Simpson and Juha-Pekka Tuovinen
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module esx_Zveg
  use CheckStops, only: CheckStop ! between point values
  use esx_Variables, only: esx, Zveg, Loc
  use Io_Routines, only : writedata

  implicit none
  private

  public  :: config_Zveg         !> choose method and params
  public  :: init_Zveg           !> allocated Zveg variables
  public  :: def_leaf_profile    !> defines LAI using beta or other functions
  public  :: Set1dPAR            !> time and z-dependent PAR
  public  :: test_Zveg           !> testing subroutine

  private ::  def_rad_prof
  private :: beta                !> beta function

  !Aug 2013: modified from EMEP land_input type, hence more complex than needed
  !for ESX, so far... 
  !> Veg type looks after the scalar data associated with one land-cover
  !!  WARNING: Much confusion with values set in the LocalVariables array

  type, public :: Veg_t
     character(len=30) :: name    = "-"
     character(len=15) :: code    = "-"
     character(len=3)  :: type    = "-"   ! Ecocystem type, see headers
     character(len=5)  :: LPJtype = "-"   ! Simplified LPJ assignment
     real    ::  hveg_max = 0.0
     real    ::  hveg     = 0.0
     real    ::  dPerh    = 0.7    ! d/hveg
     real    ::  gMax     = 0.0  ! max. conductance, m/s
     real    ::  alpha    = 0.0  !

     real :: Lm = 0.05   ! Leaf dimension (m)

     character(len=15) :: LAI_method ="input"
     character(len=15) :: profile    ="uniform"
     real :: beta_a =10.0, beta_b = 4.0  !> beta parameters, default
     real :: canopy_bottom =0.5  !> If uniform canopy, fraction of hVeg
     real :: LAI=5.0, SGS = 100, EGS = 200            ! FAKE
     real :: Astart = 100, Aend = 200
     real    ::  Eiso = 0.0   ! Emission potential isoprene, ug/g/h
     real    ::  Emtl = 0.0   ! Emission potential m-terpenes, light
     real    ::  Emtp = 0.0   ! Emission potential m-terpenes, pool
  end type Veg_t
  type(Veg_t), public  :: Veg = Veg_t()

  real, private, save :: pi        = 4*atan(1.0)
  real, private, save :: pi_div180 = 4*atan(1.0)/180.0

 contains
  !----------------------------------------------------------------------------
  !> config_Zveg chooses the method for setting LAI etc.
  !! values of the z array, hVeg, and LAI

  subroutine config_Zveg(io,writelog)
    integer, intent(in) :: io
    logical, intent(in) :: writelog 
    integer :: ilog

    namelist /esxVeg_config/ Veg
    rewind(io)
    read (io, nml=esxVeg_config)
    if(  writelog ) then
      open(newunit=ilog,file="LogConfig.Zveg")
      write(ilog,nml=esxVeg_config)
      close(ilog)
    end if
  end subroutine config_Zveg


  !----------------------------------------------------------------------------
  !> init_Zveg allocates the LAI-associated variables based upon the values of
  !! the z array, hVeg, and LAI. Should be called if veg LAI changes.

   subroutine init_Zveg()

    integer :: iz, nz, io
    real, pointer, dimension(:) :: dLAI, cumLAI

    nz = esx%nz
    dLAI => Zveg(1:nz)%dLAI
    cumLAI => Zveg(1:nz)%cumLAI

    dLAI = 0.0
    cumLAI = 0.0

    call def_leaf_profile(Veg%LAI,Veg%hVeg)

      !> Output LAI z-distribution
      open(newunit=io,file="LogLAIz.txt")
      write(io, "(a,5f7.2,i3)") "#Zveg:LAI,hveg,htrunk,a,b,nhVeg: " // &
         trim(veg%profile) , Veg%LAI,Veg%hveg,Veg%canopy_bottom, &
           Veg%beta_a, Veg%beta_b, esx%nhVeg
      write(io,"(a4,2a7,9a9)") "iz", "z", "zbnd", "dz", "dLAI", "cumLAI"
  
      do iz = esx%nhVeg, 1, -1
        write(io, "(i4,2f7.2,9f9.3)") &
         iz, esx%z(iz), esx%zbnd(iz), esx%dz(iz), dLAI(iz), cumLAI(iz)
      end do
      close(io)

   end subroutine init_Zveg

  !----------------------------------------------------------------------------
  !> Set1dveg updates vertical profiles of PAR and leaf conductance.

   subroutine Set1dPAR( )

       real :: sumPAR
       integer :: nv, nz
       logical :: daytime

       nz      = esx%nz
       nv      = esx%nhVeg
       sumPAR  = Loc%PARsun + Loc%PARshade   
       daytime = sumPAR > 1.0e-10

       Zveg(:)%PARz = sumPAR    ! above canopy (but not used!)

       Zveg(1:nv)%PARz =  def_rad_prof(sumPAR, Zveg(1:nv)%cumLAI, kRad=0.4, theta=60.0 )

       if( esx%debug_Zveg > 0 ) then
          print *, "SET1dPAR ", Zveg(1)%PARz, Zveg(nv)%PARz,Zveg(nz)%PARz
          call writedata("LogPARngs",(/"  z", "PAR", " gs"/), esx%z(1:nz), &
             reshape( (/ Zveg(1:nz)%PARz, Zveg(1:nz)%gleaf /), (/ nz, 2 /) ) &
            ,"# Vegetation radiation and condutance profiles" )
       end if

       
   end subroutine Set1dPAR

  !----------------------------------------------------------------------------

   subroutine def_leaf_profile(LAI,hVeg)
    real, intent(in)               :: LAI  &!! leaf area index (total, m2/m2)
                                         ,hVeg    !! vegetation height (m)

    real :: zPerh                  !> =z/h
    real :: fdiv                   ! hVeg*beta(a,b)
    integer :: iz, nz
    real :: a, b, zbelow, hTrunk=-999. , dzinside, sumLAI
    real, pointer, dimension(:) :: dLAI, cumLAI, z, dz, zbounds

    nz      = esx%nz
    z       => esx%z(1:nz)
    dz      => esx%dz(1:nz)
    zbounds => esx%zbnd(1:nz)
    dLAI    => Zveg(1:nz)%dLAI
    cumLAI  => Zveg(1:nz)%cumLAI

    dLAI=0.0 

    select case ( Veg%profile )
    case ( "uniform" ) 
      hTrunk = hVeg * Veg%canopy_bottom      ! 
      a = LAI/( hVeg - hTrunk )              ! LAI per unit z

      zbelow = 0.0                    
      esx%nhVeg  = -99

      ZULOOP: do  iz = 1, nz
         if( iz>1) zbelow = zbounds(iz-1)
         dzinside = 0.0
         if ( hTrunk >  zbounds(iz) ) then 
            cycle
         else if ( hVeg   <= zbounds(iz) ) then  ! canopy-top in last layer
            dzinside = hVeg-zbelow - max( hTrunk-zbelow, 0.0 )
            !print "(a,i4,8f9.2)", "Profiling h>ztop  ", iz, hTrunk, hVeg, &
            !    zbelow, zbounds(iz), hVeg-zbelow , hTrunk-zbelow, dzinside

            esx%nhVeg = iz  ! Found top layer. Will exit below

         else  
            dzinside = zbounds(iz) - max( hTrunk, zbelow ) 
            !print "(a,i4,8f9.2)", "Profiling isnide  ", iz, hTrunk, hVeg, zbelow, zbounds(iz), dzinside
         end if
         dLAI(iz) = dzinside * a / dz(iz)
         print "(a,i4,8f9.2)", "Profiling LAI: ", iz, hTrunk, hVeg, dzinside, dLAI(iz), sum(dLAI)
         if( esx%nhVeg > 0 )   exit
      end do ZULOOP
    
    case ( "beta" ) 
       a = Veg%beta_a
       b = Veg%beta_b
       fdiv = hVeg*beta(a,b)

       stop "NEEDS RE-CODE. Fails with e.g. zbnd=5,10, but hveg=6"
       ZBLOOP: do  iz = 1, nz
         zPerh = z(iz)/hVeg
         if( zPerh <= 1.0 ) then

           dLAI(iz) = LAI*zPerh**(a-1.) * (1-zPerh)**(b-1.)/fdiv

           print "(a,f6.1,i3,f6.1,3f12.4)", "beta dLAI ", hVeg, iz, z(iz), zPerh, &
                  zPerh**(a-1.) * (1-zPerh)**(b-1.)/fdiv, dLAI(iz)
         else
           print "(a,i3,2f6.1,3f12.4)", "beta xLAI ", iz, z(iz), hVeg, &
              zPerh, sum( dLAI(:)*dz(1:nz) )
            esx%nhVeg = iz  ! Found top layer. Will exit below
           exit ZBLOOP 
         end if

       end do ZBLOOP

       !esx%nhVeg = min( iz-1, nz )

     ! With beta, we need to normalise to get LAI = cumLAI(1).  

       sumLAI = dot_product ( dLAI(1:nz), dz(1:nz) ) 
       dLAI(:) = dLAI(:) * LAI/ sumLAI


    case default
       call CheckStop(.true., "esx_Zveg:unknown profile method: "//trim(Veg%profile))
    end select

   ! cumLAI is used for radiation, so we aim at centre of grid-cell
    do iz = esx%nhVeg, 1, -1
      cumLAI(iz) = dot_product ( dLAI(iz:nz), dz(iz:nz) ) - 0.5*dLAI(iz)*dz(iz)
      print "(a,i4,3f8.3)", "FinLAI ", iz, dLAI(iz), cumLAI(iz), dot_product ( dLAI(iz:nz), dz(iz:nz) )
    end do

    call CheckStop(any( cumLAI<0.0)  , "esx_Zveg:Problems cumLAI: "//trim(Veg%profile))
    call CheckStop(any(   dLAI<0.0)  , "esx_Zveg:Problems cumLAI: "//trim(Veg%profile))

   end subroutine def_leaf_profile
  
  !----------------------------------------------------------------------------
  !> Simple radiation profile in canopy
  !> @param[in]   PAR, kRad, theta
  !> @param[out]  PARz

   function def_rad_prof(I, cumLAIz, kRad, theta ) result(Iz)

    real, intent(in) :: I      &!> W/m2
                       ,kRad   &!! radiation extinction coefficient   
                       ,theta   !! angle (degrees)
    real, dimension(:), intent(in) :: cumLAIz    !> m2/m2

    real, dimension(size(cumLAIz)) :: Iz   !> W/m2

      Iz = I*exp(-kRad*cumLAIz/cos(theta*pi_div180))

   end function def_rad_prof

  !----------------------------------------------------------------------------
  !> beta function, derived using fortran gamma function. (see wikipedia)
  !  (Might move to Functions at later stage, but now just used here)

   function beta(x,y) result(b)
    real, intent(in) :: x, y
    real :: b
    b = gamma(x) * gamma(y)/gamma(x+y)
   end function beta
  !----------------------------------------------------------------------------
  !> testing code. Just call to run init_Zveg with some typical values.

  subroutine test_Zveg()
    real :: PAR = 400.0

    !call SetLeafArea(103) !> gets LAI for day 183
    call init_Zveg()

    Zveg%PARz = def_rad_prof(PAR, Zveg%cumLAI, kRad=0.4, theta=30.0 )
    !SKIP Zveg%gleaf=  def_gleaf(Zveg%PARz,gMax=0.02,alpha=0.01,nhVeg=9,daytime=.true.) ! io present triggers debug ro screen
    print *, "GLEAF ", maxval(Zveg%gsto)
    
  end subroutine test_Zveg
  !----------------------------------------------------------------------------
end module esx_Zveg

!TEST PROG AT END
!----------------------------------------------------------------------------
!SKIP! SetLeafArea gets LAI
!SKIP! Skip for now. Probably LAI should be set outisde ESX.
!SKIP
!SKIP!subroutine SetLeafArea(jday)
!SKIP!  integer, intent(in) :: jday
!SKIP!  
!SKIP!
!SKIP!    print *, "SetLeafArea ",  jday, Veg
!SKIP!    if ( Veg%LAI_method=="EMEP-latitude" ) then
!SKIP!      ! keys = max min DSGS  ESGS
!SKIP!      Veg%LAI = plf_value(reshape((/ real :: &
!SKIP!      &  Veg%SGS, Veg%LAI_vals(2),  & !min
!SKIP!      &  Veg%SGS + 30, Veg%LAI_vals(1), &
!SKIP!      &  Veg%EGS, Veg%LAI_vals(1) &
!SKIP!      &  /), (/2, 3/)), real(jday, dp))
!SKIP!    end if
!SKIP!   print *, "SET LEAF AREA ", jday, Veg%LAI
!SKIP!       
!SKIP!  end subroutine SetLeafArea
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!UNCOMMENT TO TEST JUST THIS MODULE:
!DSX program tester
!DSX   use esx_Zgrid, only : init_Zgrid, test_Zgrid
!DSX   use esx_Zveg,  only : test_Zveg, config_Zveg
!DSX   integer :: ionum
!DSX 
!DSX 
!DSX    open(newunit=ionum,file="config_esx.nml")
!DSX    call config_Zveg(ionum,writelog=.true.)
!DSX    close(ionum)
!DSX
!   If file wanted, use other ionum for TestingLAI.txt output
!
!DSX   open(newunit=ionum,file="TestingLAI.txt")
!DSX    call  test_Zgrid(ionum) 
!DSX    call  test_Zveg() 
!DSX 
!DSX   close(ionum)
!DSX end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
