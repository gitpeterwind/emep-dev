!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module EmisDef_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
implicit none
public :: EmisDef_Init    ! sets names and conv. factors of allowed emissions
public :: EmisDef_Index   ! function to find index of given pollutant name

    !----------------- basic emissions file definitions --------------------!
    !  Here we define the parameters *not* likely to change often           !
    !  between different  model versions - e.g. the size and characteristics!
    !  of emission files and sector splits.                                 !
    !-----------------------------------------------------------------------!


  !/.. First, define here the characteristics of the EMEP/CORINAIR
  !    SNAP sector emissions data-bases: 

   integer, public, parameter :: NCMAX  =  4   ! Max. No. countries per grid
!hf
   integer, public, parameter :: FNCMAX  =  3   ! Max. No. countries (with 
                                                ! flat emissions) per grid


   !/.. List of possible emissions, and their initial conversion factors:

   type, public :: emislist
      character(len=6) :: name
      real             :: conv  ! factor from emis file units to required units
   end type emislist

   integer, public, parameter :: NEMIS_DEF=8
   type(emislist), public, save, dimension(NEMIS_DEF) :: &
                       EmisDef  = emislist("notdef", -1.0)


   !/.. Sector specific information

   integer, public, parameter :: &
          NSECTORS  = 11       ! Number of SNAP-sectors in emissions
!hf
!hf SECENARIO
  integer, public, parameter :: &
          ANTROP_SECTORS=10 ! Non-natural sectors

   integer, public, parameter :: &
          ISNAP_NAT  = 11,&      ! SNAP index for volcanoe emissions 
          ISNAP_SHIP = 8         ! SNAP index for flat emissions,e.g ship
                                 ! Note that flat emissions do NOT necessarily
                                 ! belong to the same SNAP sector

!  New vertical allocation from SNAP sectors. 
!       - allocations are guesswork  - should be investigated

   integer, public, parameter :: NEMISLAYERS = 7
   real, public, parameter, &
     dimension(NEMISLAYERS,NSECTORS) ::    &

   VERTFAC =                     &  ! Vertical distribution of SNAP emissions
             reshape (           &  ! Vertical distribution of SNAP emissions
!       Ground , ...       High
    (/  0.0    , 0.00, 0.08, 0.46,0.29,0.17,0.00, & ! SNAP1= public power stations,150m nat gas
        0.5    , 0.5, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP2 = Comm./inst. combustion 
        0.0    , 0.04, 0.19, 0.41,0.30,0.06,0.0,  & ! SNAP3 = Industrial combustion !60m nat gas
        0.9    , 0.1, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP4 = Production processes
        0.9    , 0.1, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP5 = Extracton fossil fuels
        1.0    , 0.0, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP6 = Solvents
        1.0    , 0.0, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP7 = Road traffic
        1.0    , 0.0, 0.0, 0.0,0.0,0.0,0.0,       & ! SNAP8 = Other mobile (trains+planes, ...)
        !0.001  , 0.178, 0.433, 0.388,0.0,0.0,0.0,     & ! SNAP9 = Waste!60m nat gas, waste incin
        0.10  , 0.15, 0.40, 0.35,0.0,0.0,0.0,     & ! SNAP9 = Waste! .. ds + some ground level
        1.0    , 0.0, 0.0, 0.0,0.0,0.0,0.0,           & ! SNAP10= Agriculture
        1.0    , 0.0, 0.0, 0.0,0.0,0.0,0.0            & ! SNAP11= Nature
        /), &
        (/NEMISLAYERS,NSECTORS /) )!hf stakheigth

!!MADE data VERTFAC /    & ! MADE look-alike 
!!MADE       KMAX_MID  ,........      
!!MADE       Ground , ...       High
!!MADE       1.0    , 0.0,  0.0, 0.0,  & ! low sources
!!MADE       0.0    , 0.25, 0.5, 0.25    ! high sources, as MADE


contains

subroutine EmisDef_Init()

  ! Here we define the allowed emissions and their standard conversion
  ! factors. Any emissions used must be in this list, but not all  of
  ! the following are required. E.g. MADE might just use sox, nox and nh3.
  ! This will be checked in the Emissions_ml

  EmisDef(1) = emislist( "sox   ",  0.5       )   !  tonne SO2 -> tonne S 
  EmisDef(2) = emislist( "nox   ", 14.0/46.0  )   !  tonne NO2 -> tonne N
  EmisDef(3) = emislist( "co    ",  1.0       )   !  tonne CO  -> tonne CO
  EmisDef(4) = emislist( "nh3   ", 14.0/17.0  )   !  tonne NH3 -> tonne N
  EmisDef(5) = emislist( "voc   ",  1.0       )   !  tonne VOC -> tonne VOC
  EmisDef(6) = emislist( "pm25  ",  1.0       )   !  tonne pm  -> tonne pm 
  EmisDef(7) = emislist( "pm10  ",  1.0       )   !  tonne pm  -> tonne pm 
  EmisDef(8) = emislist( "pmco  ",  1.0       )   !  tonne pm  -> tonne pm 

end subroutine EmisDef_Init
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function EmisDef_Index(name) result(index)

  ! Here we find the index in EmisDef which corresponds to the pollutant name
  ! passed as an argument. If no match is found index is returned as -1.
  !---------------------------------------------------------------------------
  character(len=*), intent(in) :: name
  integer :: index, i

  index = -1
  do i = 1, NEMIS_DEF
     if ( name == EmisDef(i)%name )  then
         index = i
         exit 
     end if
  end do

end function EmisDef_Index

!_____________________________________________________________________________
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module EmisDef_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
