module esx_GetData
  use esx_Variables,  only : Loc
  use TimeDate,       only : day_of_year, date
  implicit none
  private

  public :: GetExternData

  type, public ::  InData
    type(date) :: startdate, indate
    real ::&
      tstart, timestep
    integer :: nsteps  !! Total number of data records
    integer :: nday    !! Ordinal day number
    integer :: ns      !! counter in 1.. nsteps
  end type
  type(InData), public,  save :: MetData
    
contains 
  !---------------------------------------------------------------------------
  !> Reads a test data file which contains outputs from the EMEP
  !! model. Can be used to test esx components in semi-realistic mode.
  !---------------------------------------------------------------------------

  subroutine GetExternData( metdate )

    type(date) :: metdate
    integer :: mm, dd,hh, ios, metdaynum
    integer, save :: ionum    !i/o number for input
    real :: t2C,G3_O3 !! cover,G_rh2m,L_rh, G45_O3,G3_O3,CO3ppb
    integer :: nday    !! year, ordinal day
    integer, save :: ns = 0           !! counter
    integer, save :: old_hh = -1      !! Only need to update once per hour here
    character(len=200) :: header

    if (metdate%hour == old_hh ) then
     print *, "Skip new meteo for ", metdate
     return
    end if
    print *, "Adds new meteo for ", metdate
    old_hh = metdate%hour

    Loc%latitude  = 45.0
    Loc%longitude = 15.0
    metdaynum = day_of_year( metdate%year, metdate%month, metdate%day )

    associate ( t => Loc%t2, ustar => Loc%ustar, G_rh2m => Loc%rh,  Precip => Loc%Precip, &
                 PARsun => Loc%PARsun, PARshd => Loc%PARshade )

    if( ns == 0 )then
       open(newunit=ionum,file="EsxInput_data.csv")
       read(ionum,"(a200)") header
       print *, "OPENED HEADER", trim(header)
    end if

    do while ( .true. ) 
      read(ionum,*,iostat=ios) nday,mm,dd,hh,&
           t2C,ustar,G_rh2m, PARsun,PARshd,G3_O3,Precip
           !DIRECT t2C,ustar,G_rh2m,L_rh,PARsun,PARshd,G45_O3,CO3ppb

      if ( ios /= 0 ) stop "ERROR in GetData"

      !/ Crude, but avoids need for timestamp math:
      !print *, nday, mm, dd, hh, metdaynum, metdate%month, metdate%day, metdate%hour
      if ( nday < metdaynum   ) cycle
      if ( mm < metdate%month ) cycle
      if ( dd < metdate%day   ) cycle
      if ( hh < metdate%hour  ) cycle

      ns = ns + 1
      exit
      
    end do
    print *, "GetData read PARSUN ", ns, nday,mm,dd,hh, Loc%PARsun


     MetData%nday = day_of_year(metdate%year, mm, dd) 
     MetData%indate = metdate

     t = t2C + 273.12 

    end associate
    
  end subroutine GetExternData
  
end module esx_GetData
!---------------------------------------------------------------------------!
!ESXD program test_TestData
!ESXD   use esx_GetData
!ESXD   use LocalVariables, only : L
!ESXD   integer :: n
!ESXD   real :: tstart
!ESXD 
!ESXD    do n = 1, 100
!ESXD       call GetMetData()
!ESXD       print *, n, MetData%indate, L%t2, L%PARsun
!ESXD    end do
!ESXD end program test_TestData
