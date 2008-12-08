! <CoDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
module CoDep_ml
  use CheckStop_ml, only : CheckStop
  use Par_ml,         only : MAXLIMAX, MAXLJMAX
  use Chemfields_ml, only  : so2nh3_24hr



 ! In a future version, we might test for  Ts_C > -2 instead of
 ! Ts_C > 0  as we have now.

 ! Still have RgsS_dry - in future rename to Rns_SO2_dry

  !---------------------------------------------------------------------------
  ! Calculates the  non-stomatal resistances Rns_SO2 and Rns_NH3. 
  !---------------------------------------------------------------------------
  ! For basic reference and methods, see
  !
  ! Unidoc =
  ! EMEP Report 1/2003, Part I. Unified EMEP Model Description.Valid for Rns_NH3
  !
  ! Also,
  ! RIS: Smith, R.I., Fowler, D., Sutton, M.A., Flechard, C: and Coyle, M.
  !      Atmos. Environ., 34, 3757-3777
  ! + Errata + pers. comm. with R.I.Smith
  !
  ! and
  ! Smith,  Unpublished Note
  ! Eiko Nemitz papers
  !
  !For Rns_SO2 : paper TO BE PUBLISHED SOON!
  !---------------------------------------------------------------------------

  implicit none

  public   :: make_so2nh3_24hr ! make 24hr ratio
  public   ::  CoDep_factors
  private  ::  Tabulate        ! pre-calculate many values to save CPU
 
  !/** Some parameters for the so2nh3_24hr calculations

  integer, private :: nhour 
  real, private, save ::   &             ! 24hr average ratio              
     so2nh3_hr(24,MAXLIMAX,MAXLJMAX)=1.0 !Predefined to 1.0 to make first hours
                                         !reasonable

  !/** Some parameters for the Rns calculations

  integer, private, parameter :: TMIN = -40, TMAX = 40 ! Allowed temp. range
  integer, private, parameter :: NTAB  = 100   ! No. intervals used for tabulation

   real, private, save, dimension(0:100,2) :: tab_humidity_fac
   real, private, save, dimension(0:100)   :: tab_exp_rh  ! For eqn (8.16)
   real, private, save, dimension(0:NTAB) :: &
           tab_acidity_fac,                &
           tab_F2,&                            ! For Unimod eqn (8.16)
           tab_F4, &                           ! for Rns_NH3
!hf CoDep
           tab_F3                              !For Rns SO2 

  !/ Calculated values /outputs):
   real, public, save ::  &
!hf CoDep        RgsS_dry          &!
!hf CoDep       ,RgsS_wet          &!
       humidity_fac      &! to interpolate Gns  across different RH
       ,Rns_NH3           & ! Resistance for NH3 on ground (water) surface
       ,Rns_SO2             ! Resistance for SO2 on ground (water) surface

! Resistances for SO2 in low NH3 conditions

   real, private, parameter :: CEHd = 180.0, CEHw = 100.0  !  dry, wet, m/s

   logical, private, parameter :: MY_DEBUG = .true.

contains
! =======================================================================

!hf CoDep
  subroutine CoDep_factors( so2nh3ratio24hr,so2nh3ratio, Ts_C, frh, &
                            forest, debug_flag)
! =======================================================================
!
! =======================================================================


! Input:
   
   real, intent(in) :: so2nh3ratio    ! so2/nh3 ratio
!hf CoDep
   real, intent(in) :: so2nh3ratio24hr    ! so2/nh3 ratio accumulated 
                                          ! over previous 24 hours

   real, intent(in) :: Ts_C           ! surface temp. (degrees C)
   real, intent (in):: frh             ! relative humidity (as a fraction)
   logical, intent (in):: forest      ! true if forest
   logical, intent (in):: debug_flag      ! true if forest


 ! On first call we will tabulate Rns_NH3

   logical, save :: my_first_call = .true.

  !local terms governing intermediary calculations in the evaluation of NH3_st:

   real, parameter :: BETA=1.0/22.0   ! Rns factors, see Unimod eqn (8.16)
   real :: F1, F2                     ! Rns factors, see Unimod eqn (8.16)
   real :: F3, F4                     ! Rns factors for SO2
   ! real :: F3, F4                   !  (not used)
   real :: a_SN                ! so2/nh3 after correction with 0.6
   real :: a_SN_24hr                ! so2/nh3 24hr average after correction with 0.6
!hf CoDep
   integer :: itemp            ! integer Temp in deg.C
   integer :: ia_SN            ! 10*a_SN
   integer :: ia_SN_24hr       ! 10*a_SN_24hr

   integer :: IRH              ! RH as percent
   real :: acidity_fac  ! to interpolate RgsS between high-low SO2/NH3 ratios 
 
   if ( my_first_call ) then

     call Tabulate()
     my_first_call = .false.
     write(*,*) "First CoDep call, ",  so2nh3ratio24hr,so2nh3ratio, Ts_C, frh, forest

   end if
       
   itemp     =  int( Ts_C + 0.5 )
   itemp     =  max(itemp, TMIN)   ! For safety
   itemp     =  min(itemp, TMAX)   ! For safety


!hf CoDep REVISE 0.6 factor WHY LIMIT a_SN to 2 for nh3???
    a_SN  = min(2.0,0.6*so2nh3ratio)    ! NOTE: we multiply bu 0.6 to
                              ! correct for vertical grad error in local nh3
                              ! Unidoc, eqn (8.15) 
   ia_SN = int( NTAB * a_SN/2.0 + 0.4999999 )   ! Spread values frm 0 to 2.0

!hf CoDep Cap a_SN_24hr at 3
   a_SN_24hr  = min(3.0,0.6*so2nh3ratio24hr)    ! NOTE: we multiply bu 0.6 to
                              ! correct for vertical grad error in local nh3
                              ! Unidoc, eqn (8.15) 
   ia_SN_24hr = int( NTAB * a_SN_24hr/3.0 + 0.4999999 )   ! Spread values frm 0 to 3.0


   IRH   = max( 1,  int( 100.0 * frh ) )
   if ( MY_DEBUG ) then
      if ( IRH<1 .or. IRH>100 .or. ia_SN < 0 ) then 
       print *, "CODEP ERROR ", IRH, frh, ia_SN, a_SN
       call CheckStop ( IRH<1 .or. IRH>100  , "CoDep IRH ERROR")
      end if
   end if

!hf CoDep====NOW THIS IS OLD STUFF============================
  !/ 1) Acidity factor & Rgs for sulphur:

!       acidity_fac = tab_acidity_fac( ia_SN )
!       RgsS_dry    = acidity_fac * CEHd
!       RgsS_wet    = acidity_fac * CEHw
!       call CheckStop ( RgsS_dry<0.0 .or. RgsS_wet<0.0  , "CoDep NEG ERROR")

  !/ 2) Humidity factor:  (F=forest, G=grass+other)

!       if( forest ) then
!           humidity_fac = tab_humidity_fac(IRH,1) 
!       else
!           humidity_fac = tab_humidity_fac(IRH,2) 
!       end if
!hf CoDep end=====================================================





  !/ 3) Rns_NH3   - see Unimod eqn (8.16)
  !                =RIS eq. (24), modified by CEH Note
  !     & Rns_SO2   - provisional !!!


    if (Ts_C >0 ) then    ! Use "rh" - now in fraction 0..1.0

           !F1 = 10.0 * log10(Ts_C+2.0) * exp(100.0*(1.0-frh)/7.0)
           F1 = 10.0 * log10(Ts_C+2.0) * tab_exp_rh(IRH)
           F2 = tab_F2( ia_SN  )

           Rns_NH3 = BETA * F1 * F2
           Rns_NH3 = min( 200.0, Rns_NH3)  ! After discussion with Ron
           Rns_NH3 = max(  10.0,Rns_NH3)

           !Ex F3 =  40.0 * log10(Ts_C+2.0) * exp(100.0*(1.0-frh)/25.0)
           !Ex F4 = tab_F4( ia_SN  )
           !Ex Rns_SO2 = F3 * F4

!New Formulation (article to be submitted)
! Rns_SO2_dry = 11.84  * exp(1.1*so2nh3ratio24hr) * ( frh**(-1.67) )

           F3 = tab_F3 (ia_SN_24hr) !11.84  * exp(1.1*so2nh3ratio24hr)
           F4 = tab_F4(IRH)      !frh**(-1.67)
           Rns_SO2 = F3 * F4


!Limit on So2/nh3 
!hf CoDep
          Rns_SO2 = min( 700.0, Rns_SO2)  ! Set because rh goes almost to zero occationally over the Alps, Greenland and Svalbard
                                          ! 700 chosen sort of random


           Rns_SO2 = max(  10.0,Rns_SO2) !hf CoDep SHOULD WE LIMIT IT to 10??


       if(MY_DEBUG .and. debug_flag) then
         write(*,*) "CODEP PRE rh", IRH, frh
         write(*,*) "CODEP PRE NH3", ia_SN, a_SN, F1, F2, Rns_NH3
         write(*,*) "CODEP PRE SO2",ia_SN_24hr, a_SN_24hr, F3, F4, Rns_SO2
       end if



     else if (  Ts_C > -5 ) then

!hf           Rns_NH3=200.0
           Rns_NH3 = 100.0
!hf CoDep           Rns_SO2=200.0 !Reason:Already have Rlow in Rsurface
           Rns_SO2= 100.0 ! multiply by lowTcorr in Rsurface,so effectively means 200

     else 

!hf           Rns_NH3=1000.0
           Rns_NH3= 500.0
!hf CoDep           Rns_SO2=1000.0!Reason:Already have Rlow in Rsurface
           Rns_SO2= 500.0 ! adds lowTcorr in Rsurface,so effectively means 1000. This will be for frozen surfaces

     end if !Ts_C



 end subroutine CoDep_factors

  !=======================================================================

   subroutine Tabulate()
    !/**  Tabulates humidity factors, 

     real :: a_SN, a_SN_24hr
     integer :: IRH, rh_lim, veg, ia_SN, ia_SN_24hr
     integer, parameter, dimension(2) :: &
          Rhlim = (/ 85,  75 /)   ! RH limits for F=forest, G=grass

    tab_humidity_fac(:,:) = 0.0

    ! Acidity factor

!ds bug     do ia_SN = 1, NTAB
     do ia_SN = 0, NTAB
       a_SN =  ia_SN/real(NTAB)
       tab_acidity_fac( ia_SN )  = exp( -(2.0- a_SN) )
       tab_F2 (ia_SN)            = 10.0**( (-1.1099 * a_SN)+1.6769 )
!hf CoDep       tab_F4 (ia_SN)            = 10.0**( (0.55 * a_SN)-1.0 ) 
     end do

!hf CoDep
     do ia_SN_24hr = 0, NTAB
       a_SN_24hr =  ia_SN_24hr/real(NTAB)
       tab_F3 (ia_SN_24hr)            = 11.84  * exp(1.1 * a_SN_24hr)
     enddo

     do veg = 1, 2
       rh_lim = Rhlim(veg)
       do IRH = rh_lim, 100
         tab_humidity_fac(IRH,veg) = ( (IRH-rh_lim)/(100-rh_lim) )
       end do
     end do

     do IRH = 0, 100

     enddo

     do IRH = 0, 100
          tab_exp_rh(IRH) = exp( (100.0-IRH)/7.0)
!hf CoDep  
          tab_F4(IRH)= (IRH/100.)**(-1.67)
     end do

   end subroutine Tabulate
  !=======================================================================

!=======================================================================
  subroutine make_so2nh3_24hr(hour,so2conc,nh3conc,cfac_so2,cfac_nh3)
! Calculates so2/nh3 ratio for surface, accumulated over the 
! previous 24 hours (called every hour in phychem)
! Could call more often, but then we need to keep more arrays in memory.

  implicit none

  integer, intent (in) :: hour
  integer :: i,j

  real, dimension(MAXLIMAX,MAXLJMAX), intent(in) :: so2conc
  real, dimension(MAXLIMAX,MAXLJMAX), intent(in) :: nh3conc
  real, dimension(MAXLIMAX,MAXLJMAX), intent(in) :: cfac_so2
  real, dimension(MAXLIMAX,MAXLJMAX), intent(in) :: cfac_nh3

  nhour=hour
  if (hour == 0) nhour=24
  so2nh3_hr(nhour,:,:)=so2conc*cfac_so2/max(1.0e-15,(nh3conc*cfac_nh3))
 !so2conc in mixing ratio 

  so2nh3_24hr(:,:)=0.0 !initialize each time step (hour)

  do i=1, MAXLIMAX
     do j=1, MAXLJMAX
        do nhour=1,24
           so2nh3_24hr(i,j)=so2nh3_24hr(i,j)+so2nh3_hr(nhour,i,j)/24. 
        enddo
     enddo 
  enddo

end subroutine make_so2nh3_24hr
!=======================================================================


end module CoDep_ml


