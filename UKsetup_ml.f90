module UKsetup_ml

  !/*** DESCRIPTION**********************************************************
  !/   reads in or sets data used for testing the deposition module
  !/   Most of these data are for test purposes only, for example, in respect 
  !/   of growing seasons.
  !/*************************************************************************

  use DepVariables_ml,only: NLANDUSE       &  ! No. UK land-classes
                      ,luname              &
                      ,hveg_max, b_inc, albedo, NH4_pl, SGS50, DSGS   &
                      ,EGS50, DEGS, LAImin, LAImax, SLAIlen, ELAIlen  &
                      ,g_pot_min , Sg_potlen , Eg_potlen     &
                      ,g_max     , g_min     , g_lightfac    &
                      ,g_temp_min, g_temp_opt, g_temp_max  &
                      ,RgsS      , RgsO        &
                      ,VPD_max   , VPD_min     &
                      ,SWP_max   , PWP       , rootdepth 

  use Met_ml,    only   : snow
  use Io_ml,      only   : IO_TMP, IO_HOURLY, IO_SNOW, &
                           open_file

  implicit none
  private


  public ::  ukdep_init &    ! reads in data used in UK deposition modules
           , get_growing_season

  !u7.lu    private :: ukdep_z0snow       ! reads in snow and z0 data
  !u7.lu real, public, save :: z0_nwp  ! z0 from NWP for that grid-square
  

contains

!=======================================================================
  !u7.lu subroutine ukdep_init(lat,long)
  subroutine ukdep_init(errmsg)
!=======================================================================
!   Reads in data associated with UK deposition modules, e.g. land-use
!   names, characteristics.

  !u7.lu real, intent(out) :: lat, long   ! co-ordinates from ukdep_z0_snow

  character(len=*), intent(inout) ::  errmsg
  integer :: lu      ! landuse category index
  character(len=20) ::  txt

  errmsg = "ok"  ! so far....

 ! read in biomass-associated data (LAI, etc.)
 ! **sc: NH4_pl also read in at this stage
  
      call open_file(IO_TMP,"r","ukdep_biomass.dat",needed=.true.,skip=3)

      do lu = 1, NLANDUSE
         read(unit=IO_TMP,fmt=*) luname(lu), hveg_max(lu), b_inc(lu),  &
           albedo(lu), NH4_pl(lu), & 
           SGS50(lu), DSGS(lu), EGS50(lu),DEGS(lu),  &
           LAImin(lu), LAImax(lu), SLAIlen(lu), ELAIlen(lu)
         print *, "UK biomass data ", lu, NLANDUSE, luname(lu)
      end do 
      close(unit=IO_TMP)

      albedo(:) = 0.01 * albedo(:)    ! Convert from % to fraction

   ! read in Gpot modifiers, light and temperature factors

     call open_file(IO_TMP,"r","ukdep_gfac1.dat",needed=.true.,skip=3)

     do lu = 1, NLANDUSE
       txt = luname(lu)
       read(unit=IO_TMP,fmt=*) luname(lu), g_pot_min(lu), &
          Sg_potlen(lu), Eg_potlen(lu), &
          g_max(lu), g_min(lu) , &
          g_lightfac(lu), g_temp_min(lu),g_temp_opt(lu),g_temp_max(lu)

          if ( txt /= luname(lu) ) then
             errmsg = "gfac1 problem : " // txt // luname(lu)
             return   !! Exits this loop
          endif
     end do
     close(unit=IO_TMP)

   ! read in ground surface resistance, VPD and SWP modifiers

     call open_file(IO_TMP,"r","ukdep_gfac2.dat",needed=.true.,skip=3)

     do lu = 1, NLANDUSE
       txt = luname(lu)
       read(unit=IO_TMP,fmt=*) luname(lu), RgsS(lu), RgsO(lu),  &
          VPD_max(lu), VPD_min(lu), SWP_max(lu), PWP(lu),  rootdepth(lu)

      !/ Some safety checks...

       if ( txt /=   luname(lu)    .or.  &
            PWP(lu) >=    SWP_max(lu)        ) then
               errmsg = "gfac2 problem : " // luname(lu)
               return   !! Exits this loop
       endif

     end do
     close(unit=IO_TMP)

   ! write out land-use names and numbers to help interactive start

     print *, "Available land-use classes are: "
     do lu = 1, NLANDUSE
       print "(i4,4x,a20)", lu, luname(lu)
     end do
     print *, " "

  ! read in site-specific data: snow cover, lat, long, and z0 from NWP model

  !u7.lu    call ukdep_z0snow(lat,long,z0_nwp,snow)

  end subroutine ukdep_init


!=======================================================================
    subroutine get_growing_season(lu,lat,SGS,EGS)
!=======================================================================

!   calculates the start and end of growing season for land-use
!   class "lu" and latitude "lat".  
!
!  inputs: lu, lat 
!  output: SGS, EGS
!  from module DepVariables_ml : SGS50, DSGS50

    integer, intent(in) :: lu     ! land-use index
    real, intent(in) :: lat       ! latitude (obtained via the Dep2 program
                                  ! through calling the sub-routine uk_dep_init
                                  ! from the current module) 
    integer, intent(out) :: SGS, EGS ! start and end of growing season
    !u7.lu real, intent(out) :: SGS, EGS ! start and end of growing season


      SGS = int ( 0.5 +  SGS50(lu) + DSGS(lu) * (lat-50.0) )
      EGS = int ( 0.5 +  EGS50(lu) + DEGS(lu) * (lat-50.0) )

  end subroutine get_growing_season

!=======================================================================
!u7.lu 
!u7.lu     subroutine ukdep_z0snow(lat,long,z0_nwp,snow)
!u7.lu 
!u7.lu !   reads in lat, long, NWP z0 value, and monthly snow as an index (1=present)
!u7.lu !   NB - the z0 values were originally used by the sub-grid methodology,
!u7.lu !   but aren't used anymore.  Still, just in case we read them in...
!u7.lu 
!u7.lu ! In... 
!u7.lu !    none
!u7.lu       
!u7.lu ! Out..
!u7.lu    real, intent(out) :: lat, long   ! latitude  and longitude
!u7.lu    real, intent(out) :: z0_nwp                 ! z0 from NWP
!u7.lu    integer, dimension(:), intent(out) :: snow
!u7.lu 
!u7.lu ! Local..
!u7.lu     real, dimension(7) :: class   ! from MADE-RIVM model
!u7.lu     integer :: ix, iy, imm, mm, iz0
!u7.lu 
!u7.lu      class = (/1.0e-4,1.0e-3,3.0e-1,3.0e-1,3.0e-1,3.0e-1,1.0e-3/)
!u7.lu   
!u7.lu      call open_file(IO_SNOW,"r","z0snow.in",needed=.true.)
!u7.lu 
!u7.lu      read(unit=IO_SNOW,fmt=*) ix, iy, iz0, lat, long
!u7.lu      z0_nwp = class(iz0+1)
!u7.lu       
!u7.lu       do imm = 1, 12
!u7.lu          read(unit=IO_SNOW,fmt=*) mm, snow(mm) 
!u7.lu       end do
!u7.lu   
!u7.lu       close(unit=IO_SNOW)
!u7.lu     end subroutine ukdep_z0snow
!u7.lu     !====================================================================
end module UKsetup_ml
