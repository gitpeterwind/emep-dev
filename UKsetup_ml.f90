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
                      ,f_phen_a, f_phen_b, f_phen_c, f_phen_d  &
                      ,f_phen_Slen, f_phen_Elen &
                      ,g_max     , f_min     , f_lightfac    &
                      ,f_temp_min, f_temp_opt, f_temp_max  &
                      ,RgsS      , RgsO        &
                      ,VPD_max   , VPD_min     &
                      ,SWP_max   , PWP       , rootdepth 

  use Io_ml,     only   : IO_TMP, IO_HOURLY, IO_SNOW, open_file

  implicit none
  private


  public ::  ukdep_init &    ! reads in data used in UK deposition modules
           , get_growing_season

  public :: fPhenology

contains

!=======================================================================
  subroutine ukdep_init(errmsg,me)
!=======================================================================
!   Reads in data associated with UK deposition modules, e.g. land-use
!   names, characteristics.

  character(len=*), intent(inout) ::  errmsg
  integer, intent(in) ::  me     ! processor number
  integer :: lu      ! landuse category index
  character(len=20) ::  txt

  errmsg = "ok"  ! so far....

 ! read in biomass-associated data (LAI, etc.)
 ! **sc: NH4_pl also read in at this stage
  
      !JUN06 call open_file(IO_TMP,"r","lde_biomass.dat",needed=.true.,skip=3)
      call open_file(IO_TMP,"r","JUN06_biomass.dat",needed=.true.,skip=3)

      do lu = 1, NLANDUSE
         read(unit=IO_TMP,fmt=*) luname(lu), hveg_max(lu), b_inc(lu),  &
           albedo(lu), NH4_pl(lu), & 
           SGS50(lu), DSGS(lu), EGS50(lu),DEGS(lu),  &
           LAImin(lu), LAImax(lu), SLAIlen(lu), ELAIlen(lu)
         if(me==0) write(*,*) "UK biomass data ", lu, NLANDUSE, luname(lu)
      end do 
      close(unit=IO_TMP)

      albedo(:) = 0.01 * albedo(:)    ! Convert from % to fraction

   ! read in Gpot modifiers, light and temperature factors

     call open_file(IO_TMP,"r","JUN06_gfac1.dat",needed=.true.,skip=3)

     do lu = 1, NLANDUSE
       txt = luname(lu)
       read(unit=IO_TMP,fmt=*) luname(lu), f_phen_a(lu), f_phen_b(lu), &
          f_phen_c(lu), f_phen_d(lu), f_phen_Slen(lu), f_phen_Elen(lu), &
          g_max(lu), f_min(lu) , &
          f_lightfac(lu), f_temp_min(lu),f_temp_opt(lu),f_temp_max(lu)

          if ( txt /= luname(lu) ) then
             errmsg = "gfac1 problem : " // txt // luname(lu)
             return   !! Exits this loop
          endif
     end do
     close(unit=IO_TMP)

   ! read in ground surface resistance, VPD and SWP modifiers

     call open_file(IO_TMP,"r","JUN06_gfac2.dat",needed=.true.,skip=3)

     do lu = 1, NLANDUSE
       txt = luname(lu)
       read(unit=IO_TMP,fmt=*) luname(lu), RgsS(lu), RgsO(lu),  &
          VPD_max(lu), VPD_min(lu), SWP_max(lu), PWP(lu),  rootdepth(lu)

      !/ Some safety checks...

      !These are applied for those vegetations where the VPD, SWP stuff should
      !be specified. For veg where VPD_max has been set to -1, we skip
      !these tests.

       if ( txt         /=   luname(lu)   .or.  &
            VPD_max(lu) > 0               .and.  &
                 (VPD_max(lu) >= VPD_min(lu).or. &
                  PWP(lu)     >= SWP_max(lu)     )  ) then
               errmsg = "gfac2 problem : " // luname(lu)
               return   !! Exits subroutine
       endif

     end do
     close(unit=IO_TMP)

   ! write out land-use names and numbers to help interactive start

     if( me==0) then
     print *, "Available land-use classes are: "
     do lu = 1, NLANDUSE
       print "(i4,4x,a20)", lu, luname(lu)
     end do
     print *, " "
     end if


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

    !xlat real :: xlat
    !xlat=max(30.0,lat)  ! keep functions within range of obs. used
    !xlat=min(65.0,xlat)
      


      SGS = int ( 0.5 +  SGS50(lu) + DSGS(lu) * (lat-50.0) )
      EGS = int ( 0.5 +  EGS50(lu) + DEGS(lu) * (lat-50.0) )
      EGS = max(SGS+30,EGS)  ! Safety, Just to ensure EGS > SGS!
      EGS = min(365,EGS)     ! Safety - for ca. 21N and southwards!

     ! Some limits to reflect Zhang et al's (2004) figures

      !if ( index(luname(lu),"IAM_DF") > 0 ) then
!	   SGS = max(100,SGS) 
!	   EGS = min(350,EGS) 
!      end if

  end subroutine get_growing_season

!=======================================================================
 function fPhenology(debug_flag,jday,a,b,c,d,Slen,Elen,SGS,EGS,Astart,Aend) &
 result (f)

! Input
  logical, intent(in) :: debug_flag
  integer, intent(in) :: jday
  real, intent(in) ::  a,b,c,d,Slen,Elen
  integer, intent(in):: SGS, EGS, Astart, Aend

! Output
   real :: f

        if ( jday <  SGS ) then
                f = 0.0
        else if ( jday <= Astart ) then
                f = a
        else if ( jday <= Astart+Slen ) then
                f = b + (c-b) * ( jday-Astart)/real(Slen)
        else if ( jday <= Aend-Elen ) then
                f = c
        else if ( jday < Aend ) then
                f = d + (c-d) * ( Aend-jday)/real(Elen)
        else if ( jday <= EGS ) then
                f = d
        else
                f = 0.0
        end if

        !if( debug_flag ) &
        !    write(*,"(a12,i4,2i6,f12.3)") "fPhenology",  jday, SGS, EGS, f

end function fPhenology
end module UKsetup_ml
