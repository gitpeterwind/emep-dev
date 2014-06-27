!> MODULE PointSource_ml
!!============================================================================
!! Implementation of NILU and other plume-rise methods, with Peter's
!! spread-plume
!!
!! The code will read coordinates and emissions from a file PointSources.txt.
!! Plume rise and hence effective emissions height  is calculated from the
!! stack parameters and meteorology, and the spread into the z-layers using
!! an assumed standard deviation of the emissions (=25% of the plume-rise
!! delta). Can work with any z-layers (no hard coding :-)
!!
!! Matthias Karl and Dave, Dec. 2012
!!============================================================================
!
!ToDo
!
!   Currently the meteorology assumes lowest layer has 90m depth. Should in
!   future also make wind-speed and dt/dz depend on z values.
!
!   Currently emissions make use of EmisNat arrays, since this was easier
!   for non SNAP emissions. Will modify.

module PointSource_ml

  ! --- prelim module to introduce point sources.
  !
  use CheckStop_ml,      only : CheckStop
  use ChemSpecs,         only : species
  use Functions_ml,      only : StandardAtmos_kPa_2_km, Tpot_2_T
  use GridValues_ml,     only : sigma_bnd, debug_proc &
                           , lb2ij, coord_in_processor, i_local, j_local &
                           , GridArea_m2
  use Io_Progs_ml,       only : open_file, ios, read_line, datewrite
  use LocalVariables_ml, only : Grid ! Grid-meteorology
  use MetFields_ml,      only : pzpbl, z_bnd, z_mid
  use ModelConstants_ml, only : KMAX_MID, PT, Pref, MasterProc, &
                                USES, & ! for PlumeMethod
                                DEBUG => DEBUG_EMISSTACKS
  use Par_ml,            only : me, MAXLIMAX, MAXLJMAX, limax, ljmax, &
                                IRUNBEG,JRUNBEG, & ! TMP for debug
                                gi0, gi1, gj0, gj1
  use PlumeRise_ml,      only : Plume_PreggerFriedrich, Plume_ASME, Plume_NILU
  use Setup_1dfields_ml, only : rcemis, temp, pp
  use SmallUtils_ml,     only : find_index, wordsplit
  use TimeDate_ml,       only:  nydays
  use PhysicalConstants_ml, only: AVOG,PI
  implicit none
  private

  public  :: readstacks       ! Read LPS data
  public  :: get_pointsources        ! get emissions
  private :: spread_plume     ! From Peter
  private :: gauss_integral   ! "     "
  private :: errorfunctionc   ! "     "

  integer, parameter, private :: NzMax=20   ! number of layers, tmp
  integer, parameter, private :: MAX_STACKS=9999  ! number of lines in input files
  integer, parameter, private :: MAX_POLLS =  20  ! number of emitted pollutants
  integer, parameter, private :: NeMax=10    ! number of emitted tracers
  real, private, dimension(NzMax-1), save  :: layer_z ! layer boundaries

  type, public :: stacktype
    character(len=20) :: name   ! used as label
    integer :: cc      ! country-code (could be an LPS code too in future)
    integer :: snap
    real    :: lat, long
    integer :: i,j     ! coords in EMEP domain
    real    :: hs      ! Physical ht., m
    real    :: d       ! diameter, m
    real    :: Vs      ! exit velocity
    real    :: Ts      !  stack temp (K)
    real    :: flow    ! emission
    real    :: M       ! buoyancy param
    real    :: dh4     ! Calculated plume rise at 4 m/s
    real    :: flue    ! flue gas rate (m^3/s)
    real, dimension(MAX_POLLS) :: emis
    real    :: stdev   ! std. dev in plume 
  end type stacktype
  type(stacktype), private, dimension(0:MAX_STACKS) :: stack  ! for input data

  logical, public, allocatable,dimension(:,:) :: pointsources
  integer, private, save :: nstacks = 0
  integer, private, save, dimension(MAX_POLLS) :: ispec_emis
  integer, public, save :: nemis_found

contains

!--------------------------------------------------------------------------

subroutine readstacks(io)
   integer, intent(in) :: io

   character(len=20) :: txtcode   ! Just set as CCSNAP for now
   integer :: cc, snap, i, j, iloc, jloc
   real    :: he, stdev   ! Eff. ht. and std. dev.
   integer :: n, ns, k, kk, iocode, Nz, ccmax
   real :: lat, long, x, y    !  coords
   real :: hs,Ts, diam, Vs  ! physical ht, Tempm, Diameter, flow-speed
   real :: dh4          !  Plume rise at 4m/s from Preggers+Friedrich
   character(len=100) :: fname
   character(len=200) :: txtinput  ! Big enough to contain one input record
   integer :: alloc_error
   real   , dimension(MAX_POLLS) :: emis
   integer, save, dimension(MAX_POLLS) :: icol_emis  ! column numbers of used emissions
   
   real :: he_max  ! max limit for emissions, found from max(he+2*stdev)
   real :: p
   logical :: first_header=.true., debug_flag=.false.
   character(len=20), dimension(MAX_POLLS+9) :: words
   integer :: errcode, ncol, emcol, nemis_cols, ispec, nn, nem, nwords, iemis
   debug_flag = DEBUG !! .and. MasterProc 

  ! The plume rise methodology is very approximate, so we calculate
  ! layer thickness once, for a standard atmosphere
  ! We also set k=1 for the level nearest the ground
   do kk = 2, KMAX_MID
     k = KMAX_MID - kk  +1
     p = PT + sigma_bnd(kk)*(Pref -PT)
     layer_z(k) = 1.0e3 * StandardAtmos_kPa_2_km( 0.001*p) 
     if( MasterProc) write(*,"(a,i3,2f12.3)") "StackLayer", kk, layer_z(k)
   end do

   he_max = 0.0
   ns  = 0 ! number of stack-data entry, this processor

   allocate( pointsources(MAXLIMAX, MAXLJMAX))
   pointsources(:,:) = .false.

   if ( MasterProc ) then
      fname = "PointSources.txt" !TEST
      call open_file(io,"r",fname,needed=.true.)
      call CheckStop(ios,"open_file error on " // fname )
   end if

  ! First, we read file to get number and list of country codes

   n = 0
   do 

      call read_line(io,txtinput,iocode)
      if ( iocode <  0 ) exit               ! End of file


     ! 1) Header line
      if ( txtinput(1:1) == "#" ) then  ! Figure out emissions indices
                                        ! from header
        if( .not. first_header ) cycle ! Currently units line
        first_header = .false.    

        call wordsplit(txtinput,MAX_POLLS+9,words,nwords,errcode,&
          separator=",")

        nemis_found =0     ! Total number of used emissions
        nemis_cols  = nwords - 9 ! number of emission cols in input data

        do ncol = 10,  nwords
            ispec  = find_index( words(ncol), species(:)%name ) 

            if( ispec > 0 ) then
               nemis_found = nemis_found + 1
               emcol=ncol-9  
               ispec_emis(nemis_found) = ispec
               icol_emis(nemis_found)  = emcol
               if(MasterProc) print *, "STACKem ", nemis_found, iemis, species(ispec)%name
            end if
        end do
        cycle

      else

     ! 2) Data lines
        read(unit=txtinput, fmt=*, iostat=iocode ) txtcode, cc, snap, hs, &
         lat, long, Ts, diam, Vs, ( emis(nn), nn=1,nemis_cols )
         if(MasterProc) print *, "STACK n HERE  ", n, ns, trim(txtcode),  lat, long

      end if

      n  = n + 1     ! number, this emis file
      ns = ns + 1    ! number, all emis files
      if(MasterProc) call CheckStop( n>MAX_STACKS, &
               "Too many stacks? "//trim(fname) )


     ! longitude latitude to (i,j) in any grid projection
      call lb2ij ( long, lat, x, y )

      i = nint(x)
      j = nint(y)

      iloc             = i_local( nint(x) )
      jloc             = j_local( nint(y) )
      if ((iloc>=1).and.(iloc<=limax).and.&
          (jloc>=1).and.(jloc<=ljmax)) then !on the correct processor

!TEST LATER call coord_in_processor( long, lat, iloc, jloc )

        if (debug_flag ) then
          print "(a,99i6)",'Stack Local coords are',ns,iloc,jloc
        endif

        stack(ns)%name   = trim(txtcode)
        stack(ns)%cc     = cc
        stack(ns)%hs     = hs 
        stack(ns)%d      = diam
        stack(ns)%Vs     = Vs
        stack(ns)%Ts     = Ts  
        stack(ns)%snap   = snap

        stack(ns)%lat    = lat 
        stack(ns)%long   = long

        stack(ns)%i      = iloc
        stack(ns)%j      = jloc
        pointsources(iloc,jloc)  = .true.

       ! Use air-temp = 283.15 for now, see Pregger+Friedrich
       ! Bouyancy flux
        stack(ns)%flow   =  0.25 * 9.81 *  diam**2 * Vs * (Ts-283.15)/Ts
        
       ! Flue gas flow (Pregger&Friedrich)
        stack(ns)%flue   =  0.25 * PI * Vs * diam**2

           do nn = 1, nemis_found
             ncol = icol_emis(nn)
             stack(ns)%emis(nn) = emis(ncol)
             if( MasterProc ) print *, "STACK SETS ", ns,  nn, ncol, emis(ncol)
           end do


        if(debug_flag ) write(*,"(a,i5,es10.3)") " StacksEmis ", &
           ispec_emis(1), ispec_emis(2)    !DS iemist, emisa1
        write(*,"(a,3i5,4x,2i4,4x,4i4)") " StacksCrds ", &
          me, i, j,iloc,jloc,  limax, ljmax
        if(debug_flag ) write(*,"(a,3i5,4x,2i4, 2f8.3,10f10.3)") " StacksProc ", &
         me, i, j,iloc,jloc, lat, long,  hs, dh4, he, stack(ns)%flow

        !if ( he > he_max) he_max = he
       end if

   end do 
   nstacks = ns
   if(debug_flag ) write(*,*) "Stacks Read. Done"

   if( MasterProc ) close(io)

! Find the maximum height of emissions across all  pollutants
! given as kup = k upwards from the ground
!     if ( layer_z(k) > he_max  )  then
!        Stacks%kup    = k
!        exit
!     end if
!   end do

end subroutine readstacks
!----------------------------------------------------------------------------
subroutine get_pointsources(i,j, debug_flag)
   integer, intent(in) :: i,j
   logical, intent(in), optional :: debug_flag
   integer :: n, Nz=10,k,kk,ie,iec1,iec2,iec3,iec4,iec5,iec6
   integer :: iemis,ispec, nn ! DSino,ino2,inh3
   real :: he, dh, stddev, Ta,  Ta_C, dtdz ! FAKE
   real, dimension(NzMax) :: fraction_per_layer
   real :: emiss
   real    :: Mh,ObukhovL
   real    :: uconv, uconv1
   logical :: myfirstcall = .true.

   logical :: mydebug

   mydebug = .false.
   if ( present( debug_flag )) mydebug = debug_flag

   Ta = Grid%theta_ref * Tpot_2_T( pp(KMAX_MID) )   ! 
   Ta_C = Ta - 273.15

   do n = 1, nstacks

     if (  stack(n)%i == i .and.  stack(n)%j == j ) then

       if ( USES%PlumeMethod == "ASME" ) then

          dtdz =  ( temp(KMAX_MID) - Grid%t2 ) / layer_z(1)  ! THETA or T?

          dh =  Plume_ASME( stack(n)%hs, stack(n)%flow, Grid%u_ref, Ta, dtdz )
          he = stack(n)%hs + dh
          if( mydebug .and. myfirstcall ) then
             write(*,"(a,i3,f10.4,2f10.4,a,2f10.3,es10.2)") &
            "Plume_ASME: stack he dh",n,he,dh, dtdz, &
             "Tpot?",  Grid%theta_ref, Ta, pp(KMAX_MID)
          end if

       else if ( USES%PlumeMethod == "NILU" ) then

          ObukhovL = 1.0e4  ! safety
          if (abs(Grid%invL) > 1.0e-4)  ObukhovL = 1.0/Grid%invL

          he = Plume_NILU( stack(n)%hs, Grid%ustar, stack(n)%Vs, stack(n)%d, &
             stack(n)%Ts, Ta, z_mid(i,j,KMAX_MID), pzpbl(i,j), Grid%u_ref,   &
             ObukhovL ) 

         ! For now we don't allow downwash. Too uncertain
          he = max( he, Stack(n)%hs)
          if(mydebug .and. myfirstcall ) then
            if (n==1) write(*,"(a,i3,f10.4,f10.4)") "Plume_NILU: stack he dh",  &
                  n,he,he - Stack(n)%hs
          !if (n==1) write(*,*) i,j,GridArea_m2(i,j)
          endif        

       else if ( USES%PlumeMethod == "PVDI" ) then

          ! Mh: Emitted heat flux in MW, VDI(1985)       
          Mh = 1.36e-3 * stack(n)%flue * (stack(n)%Ts-283.15)
          dh = Plume_PreggerFriedrich( Mh, Grid%u_ref, stack(n)%Vs, stack(n)%d )
          he = stack(n)%hs + dh      
          if (DEBUG .and. n==1) write(*,"(a,i3,f10.4,f10.4)") "Plume_PVDI: stack he dh",n,he,dh

       else

          call CheckStop("PlumeMethod not set! " // USES%PlumeMethod )
       end if

       dh     = he - Stack(n)%hs

       stddev = 0.25 * dh   ! Bieser et al. assumed plume top and bottom at
                            ! 0.5 dh, so we assume this was 2 std-dev
       stddev = max( stddev, 10.0 )     ! Just to force some mixing
       
       fraction_per_layer = 0.0
       call spread_plume( he, stddev, layer_z(1:Nz-1), fraction_per_layer(1:Nz),Nz)

       if( mydebug ) then
        call datewrite("STACK-"//trim(USES%PlumeMethod) //"-"//trim(stack(n)%name)&
               // " DD/MM hh:00", &  ! Format for date
         (/ he, Ta-273.15, Grid%u_ref,  Grid%invL, ObukhovL, &
              (fraction_per_layer(k), k=1,6) /), txt_pattern=.true.)
       end if

!DS       emiss(1) = stack(n)%emisnox*0.95 *uconvno      ! NO
!DS       emiss(2) = stack(n)%emisnox*0.05 *uconvno2     ! NO2
!DS       emiss(3) = stack(n)%emisnh3      *uconvnh3     ! NH3
!DS     !  emiss(4) = 31536. * uconvtr    ! tracer amine 1g/s
!DS       emiss(4) = stack(n)%emisa1 *uconvtr

       ! CONVERT EMISSION FROM kg/year to (molec/cm2/s)/dz(cm)      
       uconv=0.1 /(nydays*24.0*3600.0)                ! Kg/yr --> 10^4 g/s
       uconv=uconv/GridArea_m2(i,j)                   ! --> g/s/cm2

       do nn = 1, nemis_found
        ispec = ispec_emis(nn)
        uconv1=uconv*AVOG/species(ispec)%molwt          ! --> molecules/s/cm2
        emiss = stack(n)%emis(nn) *uconv1

        do kk=1,Nz
          k = KMAX_MID - kk + 1
          rcemis(ispec,k)  = rcemis(ispec,k) + fraction_per_layer(kk) *emiss &
                / ( 100.0*(z_bnd(i,j,k) - z_bnd(i,j,k+1)) )  ! width in cm
         if(mydebug ) write(*,"(a20,i3,f7.3,3es10.3)") " StacksEmisRC "//trim(species(ispec)%name),   &
           k, fraction_per_layer(kk), emiss, rcemis(ispec,k)
        end do !kk

       end do ! iemis 

       if(mydebug ) write(*,"(a,es10.3,es10.3,a,2es10.3)") " StacksEmisConv ",   &
         uconv,GridArea_m2(i,j),trim(species(ispec)%name),species(ispec)%molwt,emiss

     end if
   end do
   myfirstcall = .false.
end subroutine get_pointsources
!----------------------------------------------------------------------------
subroutine spread_plume( h, w, layer_z, fraction_per_layer,Nz)
!gives the fraction of a gaussian distribution between layer_z levels
!gaussian has sigma=w and is centered at h
!fraction_per_layer(1) is integral until layer_z(1)
!fraction_per_layer(Nz) is integral from layer_z(Nz-1) until infinity
  implicit none
  integer, intent(in) ::Nz !number of layers
  real, intent(in) :: h, w !plume ht, width
  real, intent(in), dimension(Nz-1) :: layer_z !(z boundaries)
  real, intent(out), dimension(Nz) :: fraction_per_layer
  real :: sum,large
  integer ::i
  
  large=1.0E10

  fraction_per_layer(1)=gauss_integral(-large,layer_z(1),h,w)

  do i=2,Nz-1
     fraction_per_layer(i)=gauss_integral(layer_z(i-1),layer_z(i),h,w)
  enddo
  fraction_per_layer(Nz)=gauss_integral(layer_z(Nz-1),large,h,w)

  !normalize (only to get last digits right)
  sum=0.0
  do i=1,Nz
     sum=sum+fraction_per_layer(i)
  enddo

  do i=1,Nz

     fraction_per_layer(i)=fraction_per_layer(i)/sum

     if( fraction_per_layer(i) > 1.0 ) then
      print *, "FRAC WRONG ", i, h,w,Nz, fraction_per_layer(i)
      call CheckStop("FRAC WRONG in spread_plume")
     end if

  enddo


end subroutine spread_plume

!--------------------------------------------------------------------------
real function gauss_integral(x1,x2,x0,sigma) result(integral)
  !integrates the area between x1 and x2 of a normalized gauss function 
  !         1/(sqrt(2PI)*sigma)*exp(-(x-x0)**2/(2sigma**2))
  
  real, intent(in)::x1,x2,x0,sigma
  real :: integ1,integ2,isigma_sqrt2

  isigma_sqrt2=1.0/(sqrt(2.0)*sigma) !scale factor for sigma
  
  call CheckStop( x2< x1, "WARNING: negative integral")
  !if(x2<x1)write(*,*)'WARNING: negative integral', x2, x1
  integ1=0.5*errorfunctionc((x1-x0)*isigma_sqrt2)
  if(x1-x0>0)integ1=1.0-integ1
  
  integ2=0.5*errorfunctionc((x2-x0)*isigma_sqrt2)
  if(x2-x0>0)integ2=1.0-integ2
  
  integral=integ2-integ1

end function gauss_integral

!--------------------------------------------------------------------------
real function errorfunctionc(x_in) result(erfc)
  !complementary error function
  !only defined for x>=0. For x<0 gives erfc(abs(x))
  !erf(0)=1 erfc(infinity)=0.0
  !gives approx 7 correct digits.

  implicit none
  real ,intent(in)::x_in
  integer ::i,j,k,n
  real ::d2x2,y,fac,PI,x

  PI=3.14159265358979323
  fac=1.0
  x=abs(x_in)

  if(abs(x)<3.0)then
     !formula for small x
     erfc=x  
     fac=-3*2
     y=1.0
     do i=1,40,1
        y=-y*x*x/i
        erfc=erfc+y*x/(2*i+1)
     enddo
     erfc=1.0-erfc*2/sqrt(pi)
!     write(*,*)'result small x',x_in,erfc
  else
     !formula for large x
     d2x2=1.0/(2*x*x)
     erfc=1.0-d2x2
     fac=-1.0
     y=d2x2
     do i=3,30,2
        fac=-fac*i
        y=y*d2x2
        erfc=erfc+fac*y
     enddo
     erfc=exp(-x*x)*erfc/(x*sqrt(pi))
!     write(*,*)'result large x',x_in,erfc
  endif
end function  errorfunctionc

!--------------------------------------------------------------------------
  

!program ReadPlume
!  use EmisDef_ml,  only :NEMIS_FILE, EMIS_FILE
!  use StackRead_ml
!  implicit none
!  integer, parameter ::Nz=20 !number of layers
!  real, dimension(Nz-1) :: layer_z !(z boundaries)
!  integer :: i
!  integer :: iemis, snap, cc
!
!  integer, parameter :: IONUM = 10
!  do i=1,Nz-1
!     layer_z(i)=50+i*50!example
!  enddo
!
!  call readstacks(IONUM, Nz, layer_z)
!
!  do iemis = 1, NEMIS_FILE
!    cc = 10
!  do i = 1, Nz-1
!    write(*,"(2a8,i4,f9.2,6f8.3)") "Final ", trim(EMIS_FILE(iemis)),  i, &
!         layer_z(i), (emis_vfac(cc,iemis,snap,i), snap=1,6)
!    write(*,"(2a8,i4,f9.2,6f8.3)") "Avg:  ", trim(EMIS_FILE(iemis)),  i, &
!         layer_z(i), (emis_vfac(99,iemis,snap,i)/(1.0e-10+sum(emis_vfac(99,iemis,snap,:))), snap=1,6)
!  end do  !i z
!  end do  !iemis 
!end program ReadPlume

end module PointSource_ml
