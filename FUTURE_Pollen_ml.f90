! <Pollen_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
module Pollen_ml

!-----------------------------------------------------------------------!
! Calculates the emission of birch pollen based on: a paper currenty 
! in prep by M. Sofiev. et al. and
! M. Sofiev et al. Towards numerical forecastin of long-range air transport of
! birch pollen: theoretical considerations and a feasibility study, 2006
! 
! Pollen emission s based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed to have 
! size aroundt 22 um in diameter and 800 kg/m3 in density. 
! 
!
!-----------------------------------------------------------------------!

!FEB2012  use ChemSpecs_tot_ml,     only : Pollen_b
  use ChemChemicals_ml,     only : species
  use EmisDef_ml,           only : QPOL,NPOL
  use Functions_ml,         only : heaviside
  use GridValues_ml ,       only : glon, glat
  use LocalVariables_ml,    only : Grid, Sub
  use MetFields_ml,         only : surface_precip, ws_10m ,rh2m
!  use CellMet_ml,           only : Grid
!  use SubMet_ml,            only : Sub, Grid !%t2 %rh
  use ModelConstants_ml,    only : KMAX_MID, KMAX_BND, &
                                   DEBUG_POLLEN,nmax, nstep,&
                                   METSTEP !, DEBUG_i,DEBUG_j
  use NetCDF_ml,            only : ReadField_CDF, Out_netCDF
  use PhysicalConstants_ml, only : AVOG ,PI
  use Par_ml,               only : MAXLIMAX,MAXLJMAX   ! => x, y dimensions
  use SmallUtils_ml,        only : find_index
  use TimeDate_ml,          only : current_date,daynumber

!-------------------------------------

  implicit none
  private

  public ::   Pollen_flux  ! subroutine

  real, public,save , allocatable,dimension(:,:,:) :: Pollen_prod
  real, public,save , allocatable, dimension(:,:)::AreaPOLL,& ! emission of pollen 
                                              heatsum,&  ! heatsum, needs to be remembered for forecast
                                              Pollen_rest! what is available pollen pr m2, fr forecast
  real, save, allocatable, dimension(:,:)   :: h, &     !Heatsum 
                                          R, &     !Summed pollen release
                                          h_day, & !Temperature summed over a day
                                          birch_frac, & ! Fraction of birch, read in
                                          h_c, &   ! Heatsum thresholds, read in
                                          corr     ! correction field for p. emission, read in
  real                                 :: scale,lim,n2m,lat_factor
  integer, save,allocatable, dimension(:,:):: day
  integer                              :: step,OpenStatus,IOSTAT,it,jt
  real,private, parameter  :: &
      T_cutoff = 273.2+3.5  & ! Cut-off temperature
     ,dH  = 50         & ! Flowering period in degree days
     ,PREC_MIN = 0.0   & ! Min cut-off precipitation [mm/h]
     ,PREC_MAX = 0.5   & ! Max cut-off precipitation [mm/h]
     ,N_TOT    = 3.7e8 & ! Could be a feld [grains/m2]
     ,RH_LOW   = 0.50  & ! Min cut--off relative humidity [%]
     ,RH_HIGH  = 0.80  & ! Min cut--off relative humidity [%]
     ,PROB_IN  = 0.2   & ! Probability for flowering to start
     ,PROB_OUT = 0.2   & ! Probability for flowering to end 
                         !(could be assumed to be larger than prob_in
     ,D_POLL   = 22.0  & ! Pollen grain diameter [um]
     ,POLL_DENS = 800.0  *1.0e3 ! Pollen density [g/m3]
   

 
   

   logical, private  :: my_first_call_pollen = .true.

  contains

!<------------------------------------------------------------------------!

   subroutine Pollen_flux (i,j,debug_flag) 
  
   implicit none

   integer, intent(in) :: i,j    ! coordinates of column
   logical, intent(in) :: debug_flag
   integer :: ipoll   !DS replaces Pollen_b
   real                :: scale,lim,invdz,n_to_mPOLL
  
! Reas in the different fields
   if ( my_first_call_pollen ) then 


      allocate(Pollen_prod(NPOL,MAXLIMAX,MAXLJMAX))
      allocate(AreaPOLL(MAXLIMAX,MAXLJMAX),heatsum(MAXLIMAX,MAXLJMAX),Pollen_rest(MAXLIMAX,MAXLJMAX))
      allocate(h(MAXLIMAX,MAXLJMAX),R(MAXLIMAX,MAXLJMAX),h_day(MAXLIMAX,MAXLJMAX))
      allocate(birch_frac(MAXLIMAX,MAXLJMAX),h_c(MAXLIMAX,MAXLJMAX),corr(MAXLIMAX,MAXLJMAX))
      allocate(day(MAXLIMAX,MAXLJMAX))

       call ReadField_CDF('pollen_data.nc','birch',birch_frac,2, &
           interpol = 'conservative',needed=.true.,debug_flag=.true.,UnDef=0.0)
       call ReadField_CDF('pollen_data.nc','h_c',h_c,2, &
           interpol = 'conservative',needed=.true.,debug_flag=.true.,UnDef=0.0)
       call ReadField_CDF('pollen_data.nc','cross',corr,2, &
           interpol = 'conservative',needed=.true.,debug_flag=.true.,UnDef=0.0)

! calculating grains to molecular weight   grains/m2*s --> mol/cm3*s 
       invdz  = 1.0e-6 / Grid%DeltaZ       ! 1/dZ [1/cm3]
       n_to_mPOLL = POLL_DENS * PI *  (D_poll**3 ) / 6.0 
       ipoll = find_index("POLLEN_B", species(:)%name ) 
       n2m = n_to_mPOLL * invdz * AVOG / species(ipoll)%molwt * 1.0e-18
       !DS n2m = n_to_mPOLL * invdz * AVOG / species(Pollen_b)%molwt * 1.0e-18
     Pollen_rest(:,:) = N_TOT
     my_first_call_pollen = .false. 
      
   end if !  my_first_call

   Pollen_prod(:,i,j) =  0.0
   AreaPOLL(i,j)      =  0.0
   if (birch_frac(i,j)==0.0 .or. h_c(i,j)==0.0) then
      scale = 0
   else 


 !  Heatsum calculations
 !  Sums up the temperatures that day for each timestep (20 min).
 !  The heatsum are then divided by timesteps, and added using the 
 !  function heatsum_calc. 
 if (nstep == 1 ) then   ! Heatsum calculation only each meteorological timestep
      if (day(i,j)==current_date%day) then
         h_day(i,j) = h_day(i,j) + Grid%t2
         step = step + 1
      else 
         day(i,j) = current_date%day
         h(i,j) = h(i,j)+ heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff)
         if (debug_flag .and. DEBUG_POLLEN) write(6,'(a10,f10.4,3f8.3,3I4)') "Heatsum: ",& 
                         h_day(i,j)/(24/METSTEP),heatsum_calc((h_day(i,j)/(24/METSTEP)),&
                         T_cutoff),h(i,j),h_c(i,j),i,j,24/METSTEP
         h_day(i,j) = Grid%t2
      end if
    ! End of heatsum calculations
end if ! heatsum calc each timestep

     heatsum(i,j) = h(i,j) ! if you want to output the heatsum field
 
      ! if heatsum is over heatsum threhold for the grid cell, the pollen
      !  emission can start calculating 
      lim = (1-prob_in)*h_c(i,j)
      
      if (h(i,j)<lim .or. Grid%t2<T_cutoff &
                     .or. surface_precip(i,j)>prec_max &
                     .or. R(i,j)>N_TOT*corr(i,j)) then 
         scale = 0

      else 

       ! reduce birch from 60 degrees north:
       if (glat(i,j) >= 60.0 ) then
          lat_factor = max(0.3, 1.0-0.005 * (glat(i,j)-60.0))
          birch_frac(i,j) = birch_frac(i,j) * lat_factor
       end if


!<----------------------------------------------------------------------------------------------- 
!  Scale is calculated from the different meteorological conditions using functions below  
         scale= corr(i,j)*(Grid%t2 - T_cutoff)/(dH*24*60*60)* &
                f_wind(ws_10m(i,j,1),Grid%wstar)* & ! wind dependece
                f_in(h(i,j),h_c(i,j),PROB_IN)* &    ! probability for flowering to start
                f_out(R(i,j),N_TOT,PROB_OUT)* &     ! probability for flowering to end
                f_cond(rh2m(i,j,1),RH_LOW,RH_HIGH)* & ! rh dependence
                f_cond(surface_precip(i,j),PREC_MIN,PREC_MAX)  ! precipitation dependence

        R(i,j)=R(i,j)+(N_TOT*scale*60*20)  ! R is to check the remains of pollen available
        Pollen_prod(QPOL,i,j) =N_tot*scale*birch_frac(i,j) !The pollen grains production
        AreaPOLL(i,j) = N_TOT*scale*birch_frac(i,j) ! [grains/m2/s]
!<-----------------------------------------------------------------------------------------------
  !!! Need to convert grains/m2*s --> mol/cm3*s
        Pollen_rest(i,j) = N_TOT - R(i,j) ! pollen grains left pr m2 after release
        Pollen_prod(QPOL,i,j) = Pollen_prod(QPOL,i,j)* n2m
   
    if (debug_flag .or. DEBUG_POLLEN ) write(6,'(a,2f7.4,f6.2)') "Result: ", &
             Pollen_prod(QPOL,i,j), R(i,j)/N_TOT,AreaPOLL(i,j)
  
      end if

    end if
  
  end subroutine Pollen_flux
  

! -------------------Functions----------------------------------------------------
  
  function heatsum_calc(t2,T_cutoff) 
   
  ! The temperature needs to be over the cutoff temperature

   real, intent(in)    :: t2,T_cutoff
   real                :: heatsum_calc  
   
   heatsum_calc = (t2-T_cutoff)*heaviside(t2-T_cutoff)

  end function heatsum_calc

!<-------------------------------------------------------------------  
  
  function f_wind(u,wstar) 

  ! Pollen emission increases with wind speeds up to 5 m/s (u_sat).
  ! This term cannot be higher than 1.5 (u_max).

  real,    intent(in)  :: u,wstar
  real                 :: u_max,u_sat,f_wind
  
  u_max = 1.5
  u_sat = 5
  
  f_wind = u_max - exp(-(u+wstar)/u_sat)
  
  end function f_wind
  
  !<--------------------------------------------------------
  
  function f_cond(x,x_min,x_max)
  
  ! This function is used for both humitidy and rain, as too much 
  ! humidity and rain stops pollen release
  
  real,    intent(in)  :: x,x_min,x_max
  real                 :: f_cond
  
  if (x>x_max) then
     f_cond = 0
  else if(x<x_min)then
     f_cond = 1
  else
     f_cond = (x_max-x)/(x_max-x_min)
  end if
  
  end function f_cond
  
  !<<------------------------------------------------------
  
  function f_in(h,h_c,prob)

  ! takes in account the uncertainity that all the trees start 
  ! to release pollen at the same time
  
  real,    intent(in)  :: H,H_c,prob
  real                 :: f_in
  
  if (h<(1-prob)*h_c) then
     f_in=0
  else if (h>(1+prob)*h_c)then
     f_in=1
  else
     f_in=(h-(1-prob)*h_c)/(2*prob*h_c)
  end if
  
  end function f_in
  
  !<<-------------------------------------------------------
  
  function f_out(R,N_tot,prob) 
  
  ! takes in account the uncertainity that all the trees stop 
  ! releasing pollen at the same time

  real,    intent(in)  :: R,N_tot,prob
  real                 :: f_out
  
  if (R<(1-prob)*N_tot) then
    f_out=1
  else if (R>(1+prob)*N_tot) then
    f_out=0
  else
    f_out = 1-(R-(1-prob)*N_tot)/(2*prob*N_tot)
  end if
  
  end function f_out
  
  
  end module Pollen_ml
  
  
  
  
  
  
  
  
  
  
  
  
  
