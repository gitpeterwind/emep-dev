module Rsurface_ml

!Changes from Version 1.0
!d1.1 - zenith angle range now to 89.9, not 85.0
!d1.3 - sinB now set to coszen, not 1/coszen

use Dates_ml,     only : daynumber, nmdays, nydays
use DepVariables_ml, only : NLANDUSE, & 
                              LU_CONIF, & ! for use of g_pot
                              g_pot, g_temp,g_vpd,g_light,g_swp,    & !u7.4
                              hveg_max  , b_inc     , albedo    ,   &
                              g_pot_min , Sg_potlen , Eg_potlen ,   &
                              g_max     , g_min     , g_lightfac,   &
                              g_temp_min, g_temp_opt, g_temp_max, &
                              RgsS      , RgsO      , RextS, RextO, &
                              VPD_max   , VPD_min   , &
                              SWP_max   , PWP       , rootdepth 
use Io_ml,        only : IO_UKDEP, open_file
!u7.lu use Metdata_ml,   only : snow, psurf, t2, METSTEP, pr 
use My_UKDep_ml, only : NDRYDEP_CALC, & ! no. of Wesely species  used
                         DRYDEP_CALC!,  & ! array with species which are needed
                         !GASNAME    !array with names of gases considered

                    
use PhysicalConstants_ml, only : KARMAN
!u7.lu - pass in met params to preserve consistency with box-model:
!u7.lu use Met_ml,   only : snow, psurf, t2, pr 
!u7.lu use Radiation_ml, only : zen, coszen, Idfuse, Idrctt,SolBio
use Radiation_ml, only : SolBio
use SoilWater_ml, only : SWP
use UKsetup_ml
use Wesely_ml, only  : Wesely_tab2, & ! Wesely Table 2 for 14 gases
                       WES_HNO3, WES_O3, & ! indices to identify HNO3, O3
                       !u7.lu Init_GasCoeff,& ! subroutine for evaluating:
                       Rb_cor,  &! correction factor used in evaluating Rb
                       DRx       !  Ratio of diffusivities to ozone
                        

use Functions_ml, only : Polygon                      
implicit none
private

public   ::  Rsurface
private  ::  g_stomatal    
private  ::  Conif_gpot
private  ::  get_glight

logical, private, save :: my_first_call = .true.
logical, private, parameter :: DEBUG_RSURF = .false.
 
    
contains
! =======================================================================

  subroutine Rsurface(lu,debug_flag,SGS,EGS, LAI,hveg,&
                      z0,ustar,Ts_C,vpd,SWP, &
                      psurf, pr, &                    !u7.lu
                      coszen, Idfuse, Idrctt, &       !u7.lu
                      snow, &                    !u7.lu
                        imm,idd,ihr,Ra_ref,g_sto,Rsur,Rb)
! =======================================================================
!
!     Description
!       calculates bulk surface resistance according to methodology
!       derived from EMEP MSC_W Note 6/00; the following pathways
!       apply for the surface resistance:
!
!        -- Rinc-- Rgs      In-canopy + soil/ground cover
!       |
!       |
!        -------- Rext      cuticular+other external surface
!       |
!       |
!        -------- Rsto      Stomatal
!
! with
!
!                           1
!  Rsur =           ____________________________
!                   LAI +  SAI  +       1
!                   ___    ___     ____________
!                   Rsto   Rext    Rinc + Rgs
!
! and we also calculate Rb here, since it (through u*) is land-use dependent 
!
!    Other gases
!
!  So far work has concentrated on ozone. However, other gases are 
!  considered using Wesely's idea of deriving resistances for ozone 
!  and SO2 first (e.g. RgsO, RgsS for Rgs) and then scaling using
!  effective Henry coefficients (H*) and reactivity coefficients (f0)
!
!
! Structure of routine
!
!  1. Calculate:
!        Rlow             low-temperature correction
!        Rinc             in-canopy resistance
!        Rsur(HNO3)  
!        Gsto(O3)         stomatal conductance (if LAI > 0)
!
!       FOR EACH remaining gas (icmp is used as an index, since cmp is assumed 
!                               to  abbreviate "component".):
!  2. Calculate ground surface resistance, Rgs
!              if (LAI<0.1) go to 4  (for snow/ice/water...)
!  3. if (LAI>0.1)  calculate Gext
!  4. Calculate Rsur(icmp)
!       END
!
! =======================================================================

!......................................
! Input:

    integer, intent(in) :: lu           ! land-use index
    logical, intent(in) :: debug_flag   ! set true if output wanted 
                                        ! for this location
!u7.lu - integer
    integer, intent(in) :: SGS, EGS     ! start and end of growing season
    real, intent(in) :: LAI             ! leaf area index (m2/m2)
    real, intent(in) :: hveg            ! height of vegetation (variable)
    real, intent(in) :: z0              ! vegetation roughness length (in m)
    real, intent(in) :: ustar           ! friction velocity (m/s)
    real, intent(in) :: Ts_C            ! surface temp. (degrees C)
    real, intent(in) :: vpd             ! vapour pressure deficit (kPa)
    real, intent(in) :: SWP             ! soil water potential (MPa)
  !u7.lu integer, intent(in), dimension(12) :: snow
    real, intent(in) ::  psurf
    real, intent(in) ::  pr
    real, intent(in) ::  coszen
    real, intent(in) ::  Idfuse
    real, intent(in) ::  Idrctt
    integer, intent(in) :: snow         ! snow=1, non-snow=0
    integer, intent(in) :: imm,idd,ihr  ! for possible use in printouts
    real, intent(in) :: Ra_ref  ! aerodynamic resistance from ca. 45 m above
                                ! the ground to (z_0+d) m above the ground

! Output:

   real, intent(out)             :: g_sto !  Stomatal conducatance (s/m)
   real,dimension(:),intent(out) :: Rsur ! bulk canopy surface resistance (s/m)
   real,dimension(:),intent(out) :: Rb   ! quasi-laminar boundary layer
                                          ! resistance (s/m)


! Working values:
   
    integer :: icmp             ! gaseous species
    integer :: iwes             ! gaseous species, Wesely tables
    real :: Hstar, f0           ! Wesely tabulated Henry's coeff.'s, reactivity
    real :: Rlow                ! adjustment for low temperatures (Wesely,
                                ! 1989, p.1296, left column) 
    real :: Rinc                ! in-canopy resistance (s/m)
    real :: xRgsS, xRgsO        ! see  DepVariables_ml
    real :: Rgs, Ggs            ! ground surface resistance + cond., any gas
    real :: Gext                ! external conductance
    real :: SAI                 ! surface area index (m2/m2)
!CORR    real :: diffc               ! molecular gas diffusivity coefficient 
!CORR                                ! (extracted from Wesely_tab2)
  ! CORR: use D_H2O, D_i instead of diffc !
    real, parameter :: D_H2O = 0.21e-4  ! Diffusivity of H2O, m2/s
                                      ! From EMEP Notes
    real            :: D_i              ! Diffusivity of gas species, m2/s

!  Just for printouts
    real :: tmpRsto, tmpRext   ! tmp,  for output (stomatal and external plant
                               ! tissue resistances, respectively)
    real :: tmpVg, tmpVg3,  tmpVgS, stoflux
    integer :: i

  ! Use SMALLSAI to decide on canopy modelling or not
    real, parameter :: SMALLSAI= 0.05  ! arbitrary value but small enough, 
                                       ! e.g. SAI usually assumes values of 
                                       ! around 3 and 4 for forests and 
                                       ! grassland, respectively
    logical :: canopy, leafy_canopy

! START OF PROGRAMME: 

  !/ Determine whether we have a canopy and if it is leafy.
  !  Surface area vegetation assumed >= LAI for tall vegetation
  !  to account for bare branches and trunk.

   SAI = LAI
   if ( hveg >=  5) then 
       SAI   = max(LAI,1.0)    ! sets min of 1 for decid.
   end if

   canopy       = ( SAI > SMALLSAI )    ! careful - can include grass
   leafy_canopy = ( LAI > SMALLSAI )    ! careful - can include grass

   if ( DEBUG_RSURF .and. debug_flag ) then
      print *, "INTO RSURF lu, sai,h", lu, SAI, hveg
   endif
   if ( DEBUG_RSURF .and. ( lu <=0 .or. lu > NLANDUSE ) ) then
      print *, "ERROR RSURF", lu, NLANDUSE, SAI, hveg, NDRYDEP_CALC
      return          
   end if

  !/** Do Rb first:
  !===========================================================================
 
  GASLOOP1: do icmp = 1, NDRYDEP_CALC
      iwes = DRYDEP_CALC(icmp)          ! index in Wesely table
                  
      if   ( hveg  >=  0.0 ) then

          Rb(icmp) = 2.0 * Rb_cor(iwes)/(KARMAN*ustar)

      if ( DEBUG_RSURF .and. Rb(icmp) < 0.0  ) then
        print *, "ERROR RSURFB", lu, icmp, iwes, Rb_cor(iwes), Rb(icmp)
        return          
      end if

      else ! water, sea, rb is calculated as per Hicks and Liss (1976)

          ! CORR - old: ! 
          ! CORR diffc = Wesely_tab2(1,iwes) ! molecular diffusivity coeff.

          ! CORR - Wesely's table gives the ratio D_H2O/D_i
          ! CORR - therefore D_i = D_H2O / Wesely

          D_i = D_H2O / Wesely_tab2(1,iwes)  ! CORR !

!CORR 19/8/2002 - spotted by pw. Corrected diffc to D_i by ds...
! should double-check equations though.
          !CORR Rb(icmp) = log( z0 * KARMAN * ustar/ diffc )
          Rb(icmp) = log( z0 * KARMAN * ustar/ D_i )
          Rb(icmp) = Rb(icmp)/(ustar*KARMAN)

         ! CORR - Rb can be very large or even negative from this
         !        formula. We impose limits:

          Rb(icmp) = min( 1000.0, Rb(icmp) )    ! CORR - - gives min 0.001 m/s!
          Rb(icmp) = max(   10.0, Rb(icmp) )    ! CORR - - gives max 0.10 m/s!

      end if
      if ( DEBUG_RSURF .and. Rb(icmp) < 0.0  ) then
        print *, "ERROR RSURFB", lu, icmp, iwes, Rb_cor(iwes), Rb(icmp)
        return          
      end if

   end do GASLOOP1
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  diffc = -99.99  ! for DEBUG
!  GASLOOP1: do icmp = 1, NDRYDEP_CALC
!      iwes = DRYDEP_CALC(icmp)          ! index in Wesely table
!                  
!      if   ( hveg  >=  0.0 ) then
!
!          Rb(icmp) = 2.0 * Rb_cor(iwes)/(KARMAN*ustar)
!
!         if ( DEBUG_RSURF .and. debug_flag ) then
!            print *, "RSURF LAND", icmp, iwes, Rb_cor(iwes), Rb(icmp)
!         end if
!      else ! water, sea, rb is calculated as per Hicks and Liss (1976)
!
!          diffc = Wesely_tab2(1,iwes)  ! molecular diffusivity coefficient
!
!          Rb(icmp) = log( z0 * KARMAN * ustar/ diffc )
!          Rb(icmp) = Rb(icmp)/(ustar*KARMAN)
!
!          if ( DEBUG_RSURF .and. debug_flag ) then
!             print *, "RSURF SEA", icmp, iwes, diffc, Rb(icmp)
!          end if
!   
!      end if
!   end do GASLOOP1


  !/** Now begin Rsur calculations:
  !===========================================================================
  !/**  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)

    Rlow = 1000.0*exp(-Ts_C - 4.0)
    Rlow = min(Rlow,9.9e4)    



!##############   1. Calculate In-Canopy Resistance, Rinc    ################
  !/** For canopies:

  !/** Calculate stomatal conductance if daytime and LAI > 0


   if( leafy_canopy  .and. &
          !u7.lu zen > 1.0e-15 .and. zen<=89.9 ) then  ! Daytime 
          coszen > 0.0 ) then  ! Daytime 

         call g_stomatal(debug_flag,lu,daynumber,imm,coszen,Idrctt,Idfuse, & !u7.lu
                         SGS,EGS, Ts_C,psurf,LAI,vpd,SWP,g_sto)
   else

         g_sto = 0.0    ! Dark, no stomatal conductance...

   end if ! leafy canopy and daytime

   if ( DEBUG_RSURF .and. debug_flag ) then
      print *, "UKDEP RSURF leafy, coszen, lai", leafy_canopy, coszen, LAI
      print "(a20,f6.1,f12.1,2i4,f12.4)", "UKDEP RSURF g_Sto", Ts_C, &
         psurf, SGS, EGS,  g_sto
   end if

  !/** Calculate Rinc, Gext 

   if(  canopy ) then   
         Rinc = b_inc(lu)* SAI * hveg  / ustar    ! Lisa's eqn. 48 for u*=0.5
         Rinc = min( Rinc, 1000.0)

   else   ! No canopy present

        Rinc = -999.9
        Gext = -999.9  !  (not needed maybe, but...)
      
   end if !  canopy


!####   2. Calculate Surface Resistance, Rsur, for HNO3 and  
!####      Ground Surface Resistance, Rgs, for the remaining 
!####      Gases of Interest                                

!In his treatment of surface resistance, Wesley does not recommend that NO_2 
! and HNO_3 be treated similarly to O_3 (see ibid. p.1296).

!The reason for deviating from Wesley here (with respect to NO_2) is that 
!whilst Wesley assumes that all NO_2 is deposited through stomata, it would 
!seem more likely that deposition on other surfaces of plants should also be 
!taken into consideration. 

 !u7.lu - previous if-test for snow replaced by more efficient
 !        multiplication, since snow=0 or 1.

  xRgsS = RgsS(lu) + Rlow  + snow * 2000.0    !u7.lu
  xRgsO = RgsO(lu) + Rlow  + snow *    0.0    !u7.lu QUERY???


!.........  Loop over all required gases   ................................

  GASLOOP2: do icmp = 1, NDRYDEP_CALC

     !-------------------------------------------------------------------------

     !  code obtained from Wesely during 1994 personal communication

        if ( DRYDEP_CALC(icmp) == WES_HNO3 ) then
            Rsur(icmp)  = max(10.0,Rlow)
            cycle GASLOOP2
        end if

     !-------------------------------------------------------------------------
     ! Calculate the Wesely variables Hstar (solubility) and f0 (reactivity)

        iwes = DRYDEP_CALC(icmp)
        Hstar =Wesely_tab2(2,iwes)    !Extract H*'s 
        f0    =Wesely_tab2(5,iwes)    !Extract f0's
    
     !-------------------------------------------------------------------------

        Rgs = 1.0/(1.0e-5*Hstar/xRgsS + f0/xRgsO)      ! Eqn. (9)
        Rgs = min(Rgs,9999.9)  
                          

     !   Use SAI to test for snow, ice, water, urban ...

   if ( DEBUG_RSURF .and. debug_flag ) then
      print "(a20,2i3,es12.3,4f9.3)", "UKDEP RSURF LOOP2", icmp, iwes, &
         Hstar, f0, xRgsS, xRgsO, Rgs
   end if
       if ( canopy  ) then   

         ! ###   3. Calculate Cuticle conductance, Gext   ################
         ! ###      and  Ground surface conductance Ggs:

         ! Corrected for other species
         ! using Wesely's eqn. 7 approach. (We identify leaf surface resistance
         ! with Rext/SAI.)

           Gext  = 1.0e-5*Hstar/RextS + f0/RextO
  
           Ggs = 1.0/( Rgs + Rinc )


         ! ##############   4. Calculate Rsur for canopies   ###############

           Rsur(icmp) = 1.0/( LAI*DRx(iwes) *g_sto + SAI*Gext + Ggs )
           if ( Rsur(icmp) < 1.0 ) then
             print *, "LOWRSUR", icmp, iwes, Rsur(icmp), lu
             print *, "LOWRSUR bits", LAI*DRx(iwes) *g_sto, SAI*Gext, Ggs
             print *, "LOWRSUR bit2", LAI,DRx(iwes),g_sto, SAI
             print *, "LOWRSUR H,f0", Hstar,f0
             print *, "LOWRSUR Rx",  RextS, RextO    
             print *, "LOWRSUR h,ustar,cosz",  hveg, ustar, coszen
             print *, "LOWRSUR gs",  Rgs, Rinc, RgsS(lu), Rlow, xRgsS
           end if
              
   if ( DEBUG_RSURF .and. debug_flag ) then
      print "(a20,2i3,es10.3,4f9.3)", "UKDEP RSURF CANOPY", icmp, iwes, &
           Hstar, f0, Gext, Ggs, Rsur(icmp)
   end if
       

       else   ! Non-Canopy modelling:

           Rsur(icmp) = Rgs
           Ggs = 1.0/ Rgs
   if ( DEBUG_RSURF .and. debug_flag ) then
      print "(a20,2i3,3f12.4)", "UKDEP RSURF NON-C", icmp, iwes, Rsur(icmp), Hstar, f0
   end if

        end if  ! end of canopy tests 

    !  write out resistances and Gsur  for comparison.
     
        if( DRYDEP_CALC(icmp) == WES_O3 ) then
          
           if (g_sto > 0.0) then 
              tmpRsto = 1.0/(LAI*DRx(iwes)*g_sto)
           else    !ds1   if (g_sto <= 0.0) then
              tmpRsto = -99.0
           end if

           if( canopy )  then 
              tmpRext = 1.0/(SAI*Gext)
           else 
              tmpRext = -99.9
           end if
   
           if (my_first_call ) then  !  Header


             if ( DEBUG_RSURF .and. debug_flag ) then
                 call open_file(IO_UKDEP,"w","Rsurf.out",needed=.true.)
                 write(unit=IO_UKDEP, fmt="(3a3,a5,a6,7a7,a9)") &
                 "mm", "dd", "hh", &
                 "Ts","u*","Ra","Rb", "Rgs", "Rinc", "Rext", "Rsto", &
                 "Rsur","  (=>cm/s)"
             end if
             my_first_call = .false.
           end if ! first_call 

             if ( DEBUG_RSURF .and. debug_flag ) then
                print *, "DEBUG_RSURF LAI, SAI, hveg " , LAI, SAI, hveg
                print "(3i3,f5.0,f6.2,2f7.1,5f7.0,f9.3)",  &
                 imm, idd, ihr, Ts_C, ustar, &
                Ra_ref, Rb(icmp), Rgs, Rinc, tmpRext, tmpRsto, &
                Rsur(icmp), 100.0/(Ra_ref + Rb(icmp) + Rsur(icmp))
                write(unit=IO_UKDEP,fmt="(3i3,f5.0,f6.2,2f7.1,5f7.0,f9.3)")  &
                 imm, idd, ihr, Ts_C, ustar, &
                Ra_ref, Rb(icmp), Rgs, Rinc, tmpRext, tmpRsto, &
                Rsur(icmp), 100.0/(Ra_ref + Rb(icmp) + Rsur(icmp))
            end if
        end if      ! WES_O3
   
   
           !-------------------------
           !call monthly_collect(ihr,LAI,DRx,g_sto,SAI,Gext,Ggs, &
           !          Ra_ref,Rb,Rsur(icmp))     !This call MIGHT be used.
           !-------------------------
       
      !ds1  end if

  end do GASLOOP2

 end subroutine Rsurface

!=======================================================================

    subroutine g_stomatal(debug_flag,lu,jday, &
                         imm,coszen,Idrctt,Idfuse, & !u7.lu
                         SGS,EGS,Ts_C,pres,LAI,vpd,SWP,g_sto)
!=======================================================================

!    Calculates stomatal conductance g_sto based upon methodology from 
!    EMEP MSC-W Note 6/00:
!
!    Gsto = [g_max * g_pot * g_light * g_temp * g_vpd * g_smd ]/41000.0
!
! Inputs:

    logical, intent(in) :: debug_flag   ! set true if output wanted 
  integer, intent(in) :: lu           ! land-use index (max = nlu)
  integer, intent(in) :: jday         ! days since 1st Jan.
  integer, intent(in) :: imm          ! month - (just for output)
  real, intent(in) :: coszen          ! cos of zenith angle
  real, intent(in) ::  Idfuse         ! Diffuse radn, u7.lu
  real, intent(in) ::  Idrctt         ! Direct  radn, u7.lu
  integer, intent(in) :: SGS, EGS     ! start, end of growing season (day no.)
  real, intent(in) :: Ts_C            ! surface temperature at height 
                                      ! 2m (deg. C)
  real, intent(in) :: pres            ! surface  pressure (Pa) 
  real, intent(in) :: LAI             ! leaf area index   (m2/m2)
  real, intent(in) :: vpd             ! vapour pressure deficit (kPa)
  real, intent(in) :: SWP             ! soil water potential (MPa)

! Outputs:
 
  real, intent(out) :: g_sto         ! stomatal conductance

!u7.4 real :: g_pot                         ! stomatal conductance age factor 
!u7.4 real :: g_temp                        ! stomatal conductance temperature factor
!u7.4 real :: g_vpd                         ! stomatal conductance vpd factor
!u7.4 real :: g_light                       ! stomatal conductance light factor
!u7.4 real :: g_swp                         ! stomatal conductance soil water factor

        
!..1 ) Calculate g_pot. Max value is 1.0.
!---------------------------------------
!u7.lu - make SGS, EGS real here, to keep Polygon function general

        g_pot =  Polygon(jday,nydays,g_pot_min(lu),g_pot_min(lu),1.0,&
                         real(SGS), Sg_potlen(lu), real(EGS), Eg_potlen(lu))
!la
!la     For coniferous forest we need to use a different type of 
!la     growing season to get g_pot, with SGS=120, EGS=485.
!la:
      if ( lu == LU_CONIF )  then
           g_pot =  Polygon(jday,nydays,g_pot_min(lu),g_pot_min(lu),1.0, &
                         120.0, Sg_potlen(lu), 485.0, Eg_potlen(lu)) 
           call Conif_gpot(imm,g_pot) !The subroutine conif_g_pot 
                                      !is defined below       
      end if

      ! Note that g_pot_min=g_pot at SGS.
  


!..2 ) Calculate g_light 
!---------------------------------------

  call get_glight(coszen,Idrctt,Idfuse,g_lightfac(lu),LAI,albedo(lu),g_light)    
  
  !The subroutine get_glight is defined below.

!..3) Calculate  g_temp
!---------------------------------------

!ds1  g_temp = 2.0
  
  g_temp = (Ts_C - g_temp_opt(lu)) / &
           (g_temp_opt(lu) - g_temp_min(lu))
  g_temp = 1.0 - (g_temp*g_temp)
  !u7.lu g_temp = max(g_temp,g_min(lu) )


!..4) Calculate g_vpd
!---------------------------------------

 g_vpd = g_min(lu) + (1.0-g_min(lu)) * (VPD_min(lu)-vpd )/ &
                                    (VPD_min(lu)-VPD_max(lu) )
 g_vpd = min(g_vpd, 1.0)
 !u7.lu g_vpd = max(g_vpd, g_min(lu))


!..5) Calculate g_swp
!---------------------------------------

  !/  Use SWP_Mpa to get g_swp. We just need this updated
  !   once per day, but for simplicity we do it every time-step.

       g_swp = g_min(lu)+(1-g_min(lu))*(PWP(lu)-SWP)/(PWP(lu)-SWP_max(lu))
       g_swp = min(1.0,g_swp)
       !u7.lu g_swp = max(g_min(lu),g_swp)


!.. And finally,
!---------------------------------------
! (using factor 41000 from mmol O3/m2/s to s/m given in Jones, App. 3
!  for 20 deg.C )

   g_sto = (g_light * g_temp * g_vpd * g_swp )
    !if(  debug_flag  .and. lu == 10 ) then
    !   print "(a10,5f12.4)", "MOOR PRE", g_temp, g_light, g_vpd, g_swp, g_sto
    !end if 
   g_sto = max( g_sto,g_min(lu) )
   g_sto = ( g_max(lu) * g_pot * g_sto )/41000.0
    !if(  debug_flag  .and. lu == 10 ) then
    !   print "(a10,3f12.4)", "MOOR Post", g_sto, g_max(lu)/41000.0, g_pot
    !end if 

   if ( g_sto < 0 ) then
     print *, "GSTO NEG" , jday, g_sto, g_pot, g_light, g_temp, g_vpd
   end if


  end subroutine g_stomatal


!la =====================================================================
    subroutine Conif_gpot(imm,g_pot)
!la =====================================================================
!la   modifies g_pot (g_age) for effect of older needles, with the simple
!la   assumption that g_age(old) = 0.5.
!
   !/ arguments

    integer, intent(in) :: imm    ! month
    real,   intent(inout) :: g_pot   ! Requires initial input of g_pot 
                                     ! (once obtained as output from g_stomatal)

   !/ Some parameters:
   !  Proportion of needles which are from current year:
    real, parameter, dimension(12) :: Pc = (/  &
                   0.53, 0.535, 0.54, 0.545, 0.1, 0.15,  &
                   0.27,  0.36, 0.42,  0.48, 0.5,  0.5  /)

    real, parameter :: G_POTOLD = 0.5  ! value of g_pot for old needles



!needles from current year assumed to have g_pot as evaluated above;
!needles from previous years assumed to have g_pot of 0.5
!The sum of the g_pot's for the current year is added to the sum of the
!g_pot's for previous years to obtain the overall g_pot for the landuse 
!category temp. conif. forests.

    g_pot = Pc(imm)*g_pot + (1.0-Pc(imm))*G_POTOLD

  end subroutine Conif_gpot
! =====================================================================

!===========================================================================
    subroutine get_glight(coszen,Idrctt,Idfuse,g_lightfac,LAI,albedo,g_light)
!===========================================================================
!
!    Calculates g_light, using methodology as described in Emberson et 
!    al. (1998), eqns. 31-35, based upon sun/shade method of  
!    Norman (1979,1982)

!     input arguments:

    real, intent(in) :: coszen    ! cos(zen), zen=zenith angle
    real, intent(in) :: Idrctt    ! total direct solar radiation (W/m^2)
    real, intent(in) :: Idfuse    ! diffuse solar radiation (W/m^2)
    real, intent(in) :: g_lightfac ! land-use specific light factor
    real, intent(in) :: LAI       ! leaf area index (m^2/m^2), vegetation        
                                  ! specific
    
    real, intent(in) :: albedo    ! Fraction, 0-1

!     output arguments
  
    real, intent(out) :: g_light    

!   Some parameters

    real, parameter :: &
      PARfrac = 0.45,   & ! approximation to fraction (0.45 to 0.5) of total 
                          ! radiation in PAR waveband (400-700nm)
      Wm2_uE  = 4.57,   &  ! converts from W/m^2 to umol/m^2/s
      cosA    = 0.5       ! A = mean leaf inclination (60 deg.), where it is 
                          ! assumed that leaf inclination has
                          ! a spherical distribution

!     internal variables:

    
    real :: sinB      ! B = solar elevation angle  = complement of zenith angle
    
    real :: LAIsun    ! sun-fraction of LAI
    real :: LAIshade  ! shade-fraction of LAI
    real :: PARsun    ! sun-fraction of PAR
    real :: PARshade  ! shade-fraction of PAR
       
    real :: g_sun    ! sun-fraction contribution to g_light
    real :: g_shade  ! shade-fraction contribution to g_light


    !DEP1.2 error sinB = 1.0/coszen  ! = cos(zen)
    sinB = coszen  ! = cos(zen)


    LAIsun = (1.0 - exp(-0.5*LAI/sinB) ) * sinB/cosA
    
    LAIshade = LAI - LAIsun

! PAR flux densities evaluated using method of
! Norman (1982, p.79): 
! "conceptually, 0.07 represents a scattering coefficient"  

    PARshade = Idfuse * exp(-0.5*LAI**0.7) +  &
               0.07 * Idrctt  * (1.1-0.1*LAI)*exp(-sinB)   

    PARsun = Idrctt *cosA/sinB + PARshade

!.. Convert units, and to PAR fraction
!.. and multiply by albedo

    PARshade = PARshade * PARfrac * Wm2_uE * ( 1.0 - albedo )
    PARsun   = PARsun   * PARfrac * Wm2_uE * ( 1.0 - albedo )


    g_sun   = (1.0 - exp (-g_lightfac*PARsun  ) ) * (LAIsun  /LAI)
    g_shade = (1.0 - exp (-g_lightfac*PARshade) ) * (LAIshade/LAI)

    g_light = g_sun + g_shade

  end subroutine get_glight

!--------------------------------------------------------------------

end module Rsurface_ml
