!> <AOTx_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **********************************************************************! 

module AOTx_ml
  use CheckStop_ml,  only : checkStop, StopAll
  use Chemfields_ml, only : xn_adv, cfac
  use ChemSpecs,     only : IXADV_O3
  use DO3SE_ml
  use GridValues_ml, only : debug_li, debug_lj, debug_proc, i_fdom, j_fdom
  use Io_Progs_ml,   only : datewrite
  use LandDefs_ml,   only : LandType
  use LocalVariables_ml, only : Grid, L
  use MetFields_ml, only: zen
  use ModelConstants_ml, only : dt_advec, KMAX_MID, DEBUG  & 
     ,PPBINV ! 1.0e9, for conversion from mixing ratio to ppb
  use NumberConstants, only : UNDEF_R, UNDEF_I
  use OwnDataTypes_ml, only : TXTLEN_DERIV, TXTLEN_SHORT
  use Par_ml, only : MAXLIMAX, MAXLJMAX, limax, ljmax, me
  use TimeDate_ml, only : current_date, print_date, jday => effectivdaynumber
  implicit none
  private

  public :: Calc_AOTx          ! called from AddMosaicOutput, My_DryDep
  public :: Calc_POD           ! called from AddMosaicOutput
  public :: Calc_GridAOTx      ! called from Derived_ml
  public :: Calc_SPOD          ! Experimental, JAN2013

! Limit of daylight zenith angle for AOTs
  integer, private,  parameter :: AOT_HORIZON  = 89         
  integer, public, parameter:: STARTMONTH_FOREST=4,ENDMONTH_FOREST=9&
                                ,STARTMONTH_CROPS=5,ENDMONTH_CROPS=7 ! EU only!
  logical, parameter, private :: T=.true., F=.false.

! VEGO3 definitions for PODY(formerly AFstY) and AOTX
!
! PODY vs Fst - not that the accumulated PODY  is processed here.
! the instantaneous Fst is set as for Canopy O3 in METCONCS
          ! N.B. AOTs have several definitions. We usually want
          ! the ICP-veg Mapping Manual (MM) ones. Other
          ! possibilities are EU (8-20daytime) or UN (May-July for
          ! crops)
   !================== 

    type, public:: O3cl_t
       character(len=TXTLEN_DERIV) :: name = '-' ! e.g. POD1_IAM_DF
       character(len=TXTLEN_SHORT) :: class = '-' ! POD or AOT
       real    :: Threshold = UNDEF_R     ! Threshold or CL, e.f. AOTx or AFstY
       character(len=TXTLEN_SHORT) :: defn = '-' !  MM or EU definitions
       character(len=TXTLEN_SHORT) :: txtLC = '-' !  CF, DF, IAM_CF etc.
       logical :: RelSGS = .false.   ! true if accumulation period is relative to
                              ! start of growing season (SGS)
                              ! can be false it fixed, e.g. April 1st
       integer :: SAccPeriod = UNDEF_I  ! Start of accumulation period, either rel 
                              ! to SGS or day number, days
       integer :: EAccPeriod = UNDEF_I  ! End ...., days
       integer :: iotype = UNDEF_I  ! 2=>IOU_YEAR, 3=>IOU_MON, 4=>IOU_DAY, 5-7=>IOU_HOUR
    end type 

   integer, parameter, private :: MAX_NUM_VEGO3 = 60
   type(O3cl_t), public, save, dimension(MAX_NUM_VEGO3) :: OutputVegO3 = O3cl_t()
   integer, save, public :: nOutputVegO3 = 0

    type(O3cl_t), public, allocatable, dimension(:) :: &
     VEGO3_OUTPUTS

 ! Grid%surf_o3_ppb = o3_45 * cfac, before deploss
 ! Grid%surf_o3_ppb1 = o3_45 * cfac, after deploss
 ! cano3_ppb from O3                before deploss

contains
 !=========================================================================
 ! Calc_AOTx called from MosaicOutputs, at end of DryDep calculation.
 ! after deploss
  subroutine Calc_AOTx(iO3cl,iLC, aot )
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: aot

    real    :: o3
    integer :: i,j, mm,hh
    real :: X
    logical :: dbg, is_MM, is_EU, is_XX
    logical, save :: first_call = .true.
    character(len=*),parameter :: sub='CalcAOTx:'
    character(len=50) :: txt

   ! MM (Mapping Manual) means use daylight O3, EU uses 7:00 -- 18:59 UTC

    is_EU =  VEGO3_OUTPUTS(iO3cl)%defn == 'EU' .or. VEGO3_OUTPUTS(iO3cl)%defn == 'XEU'
    is_MM =  VEGO3_OUTPUTS(iO3cl)%defn(1:2) == 'MM' .or. &  ! 1:2 allows MM:Plus5
             VEGO3_OUTPUTS(iO3cl)%defn == 'XMM'
    is_XX =  VEGO3_OUTPUTS(iO3cl)%defn == 'XMM' .or. VEGO3_OUTPUTS(iO3cl)%defn == 'XEU'
    mm = current_date%month
    hh = current_date%hour

    i = Grid%i
    j = Grid%j
    X = VEGO3_OUTPUTS(iO3cl)%Threshold
    dbg =  DEBUG%AOT .and. debug_proc .and. i == debug_li .and. j == debug_lj

    if( dbg ) then
      txt = sub// trim(VEGO3_OUTPUTS(iO3cl)%name)// print_date()
    end if

    aot = 0.0

    if (  is_XX .and. LandType(iLC)%is_forest .and. &
         ( mm<STARTMONTH_FOREST .or. mm>ENDMONTH_FOREST)   )  then  
       ! for both EU and MM-EMEP so far. Reconsider for MM
          if( dbg) print *, txt//" RETURN jd, SGS-EGS ", jday, L%SGS,L%EGS
          RETURN
    end if

    if( LandType(iLC)%is_forest .and. dbg ) then
          write(*,*) txt//"DEBUGOT CONT jd, SGS-EGS,smm ", &
             jday, L%SGS,L%EGS, STARTMONTH_FOREST
    end if

   ! If night, or outside growing season, we simply exit with aot=0
   ! EU AOT for 8:00 -- 20:00 CET, is really 8:00 -- 19:59 CET 
   ! Or:        7:00 -- 18:59 UTC
   ! (nb hour is integer value)

    if ( is_EU  .and. (hh  <  7 .or. hh > 18 )) RETURN
    if ( is_MM .and.  Grid%izen >= AOT_HORIZON ) RETURN !UN or MM use daylight

  ! the wheat growing season is based upon the latitude function
    if ( is_MM ) then
      if ( dbg ) write(*,*) txt//'AOTiO3CL ', iLC, iO3cl, L%SGS, jday,  &
           LandType(iLC)%is_crop, LandType(iLC)%is_iam

      if ( jday < L%SGS .or. jday > L%EGS ) then
         RETURN
      end if

      o3 = L%cano3_ppb   ! canopy top O3 for MM

    else if ( is_EU ) then
      if (   LandType(iLC)%is_crop .and. &
            ( mm<STARTMONTH_CROPS .or.mm>ENDMONTH_CROPS) ) then
           RETURN
      end if

      o3 = Grid%surf_o3_ppb  ! Use 3m O3 for EU

      else 
        if ( first_call ) call CheckStop(txt//"AOT not MM or EU!")
    end if

   ! ======== Plus and Minus 5 ppb used for sens. tests
      if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Plus5' ) then  ! nmole/m2/s
         o3 = o3 + 5
      else if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Minus5' ) then  ! nmole/m2/s
         o3 = max( 0.0, o3 - 5 )
      end if
    !========== Calculate ========================

    if ( o3>X ) then  ! Add AOT, scaling for time-fraction

        aot = (o3-X)  ! dt_advec takes care of 3600s elsewhere

    end if

    if (  dbg ) call datewrite(sub//"AOTxdebug" // "defn:" // &
      trim(VEGO3_OUTPUTS(iO3cl)%defn), iLC,  &
      (/ real(Grid%izen), X,  o3, aot /) )

    first_call = .false.

  end subroutine Calc_AOTx

 !=========================================================================
 ! Calc_GridAOTx called from Derived, after deploss

  function Calc_GridAOTx( iX ) result (aot)

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    integer, intent(in) :: iX  ! usually 40 (ppb)

    real, dimension(limax,ljmax) :: aot

    integer :: izen                    ! integer of zenith angle
    real :: o3, o3_ref                 ! Ozone (ppb) - needed if AOTs
    integer :: i, j 
    logical :: dbg
    character(len=*), parameter :: dtxt='CalcGridAOT:'


    aot = 0.0

    !--------- ONLY April-Sept -------------------------
!GINA     if(   current_date%month<STARTMONTH_FOREST&
!GINA       .or.current_date%month>ENDMONTH_FOREST) return
    !--------- ONLY April-Sept -------------------------
!Sep 13th, keep for now. Irrelevant for global studies though
    if(  current_date%month<STARTMONTH_FOREST&
     .or.current_date%month>ENDMONTH_FOREST   )  RETURN

    do i=1, limax
        do j=1, ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then                 ! Use 3m O3 for EU

                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV

                aot(i,j) = max( o3 - iX , 0.0 )   ! Definition of AOTs

           end if
        end do
    end do
    if ( DEBUG%AOT .and. debug_proc ) then
       i=debug_li; j=debug_lj
       o3_ref = xn_adv(IXADV_O3,i,j,KMAX_MID) * PPBINV
       o3     = o3_ref * cfac(IXADV_O3,i,j)
       call datewrite(dtxt, (/ zen(i,j), o3_ref, o3, aot(i,j) /) )
    end if

  end function Calc_GridAOTx


!=========================================================================
!Calc_POD called from DryDep -> Add_MosaicOutput,  after deploss

  subroutine Calc_POD(iO3cl,iLC, pod, debug_flag, debug_txt )
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: pod

    real    :: o3, o3_3m, o3_3mC, pod0
    integer :: i,j, spod, epod
    logical :: dbg
    real :: Y

    pod = 0.0
    dbg = ( DEBUG%AOT .and. debug_flag )

    if ( Grid%izen >= AOT_HORIZON ) then  !UN or MM use daylight
        if ( dbg ) write(*,*) "PODxZen ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, Grid%izen, L%SGS, L%EGS
        return
    end if

    ! Start and end of POD accumulation period:
    ! ( Some PODs are only for a specific period, different to 
    !   growing season (SGS), notably TC (wheat) which has 55 days.)

    if ( VEGO3_OUTPUTS(iO3cl)%RelSGS .eqv. .true. )  then
      spod = L%SGS + VEGO3_OUTPUTS(iO3cl)%SAccPeriod
      epod = L%SGS + VEGO3_OUTPUTS(iO3cl)%EAccPeriod
    else
      spod = L%SGS
      epod = L%EGS
    end if
    if ( dbg ) write(*,'(2a,5i5)') "PODuZen ",&
     trim(VEGO3_OUTPUTS(iO3cl)%name)//'def:'//trim(VEGO3_OUTPUTS(iO3cl)%defn),&
      iO3cl, jday, Grid%izen, spod, epod
    if ( jday < spod .or. jday > epod ) then
       if ( dbg ) write(*,*) "PODxJday RETURN "
       return
    end if

    i = Grid%i
    j = Grid%j
    Y = VEGO3_OUTPUTS(iO3cl)%Threshold ! nmole/m2/s

    o3 = L%cano3_ppb   ! canopy top O3 for MM
    o3_3m = Grid%surf_o3_ppb 
    o3_3mC = 0.5*(  Grid%surf_o3_ppb  + Grid%surf_o3_ppb1 ) ! half-way in dep
    if( dbg) print *, "PODCOMP ", o3_3m, xn_adv(IXADV_O3,i,j,KMAX_MID) * PPBINV * cfac(IXADV_O3,i,j), &
      Grid%surf_o3_ppb, Grid%surf_o3_ppb1  ! Use 3m O3 for EU

   ! Add fluxes if Y exceeded:

     pod  = max(L%FstO3 - Y,0.0)

  !GINA TEST
    pod0 = pod
    !! o3_3m = -999.0
    if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Plus5' ) then  ! nmole/m2/s
       if( o3_3m>0.1 ) pod = pod * (o3_3m+5)/o3_3m
   if ( dbg ) write(*,"(a,i4,5f12.5)") "PODaaaa1 "//&
           trim(VEGO3_OUTPUTS(iO3cl)%name)//':'// trim(VEGO3_OUTPUTS(iO3cl)%defn),&
           iO3cl, o3_3m, L%FstO3, (o3_3m+5)/o3_3m, pod, pod0
 
    else if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Minus5' ) then  ! nmole/m2/s
       o3_3m = xn_adv(IXADV_O3,i,j,KMAX_MID) * PPBINV * cfac(IXADV_O3,i,j)
       if( o3_3m>0.1 ) pod = pod * max(o3_3m-5, 0.0)/o3_3m 
   if ( dbg ) write(*,*) "PODaaaa2 ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, trim(VEGO3_OUTPUTS(iO3cl)%defn)
    else if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:C' ) then  ! nmole/m2/s
       if( o3_3m>0.1 ) pod = pod *  o3_3mC/o3_3m
    else
        if ( dbg ) write(*,*) "PODaaaa3 ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, trim(VEGO3_OUTPUTS(iO3cl)%defn)
    end if
   !END GINA TEST
       

    if ( dbg ) then
       write(6,"(a,L2,3i4)") "PODACC ", VEGO3_OUTPUTS(iO3cl)%RelSGS, &
           VEGO3_OUTPUTS(iO3cl)%SAccPeriod, VEGO3_OUTPUTS(iO3cl)%EAccPeriod, &
           jday
       write(6,"(a,L2,3i4)") "PODPM  ", VEGO3_OUTPUTS(iO3cl)%RelSGS, &
           VEGO3_OUTPUTS(iO3cl)%SAccPeriod, VEGO3_OUTPUTS(iO3cl)%EAccPeriod, &
           jday
       call datewrite("YYY-debug_pod-defn:" // &
         trim(VEGO3_OUTPUTS(iO3cl)%defn), iLC, &
           (/ real(Grid%izen), Y,  o3, L%FstO3, pod /) )
       if(current_date%hour > 10 .and. pod > 0.1 ) then
       write(6,"(a,2f9.4,2es12.4)") "PODCOMP ",  o3, o3_3m, pod0, pod
            !call StopAll('TMP DSEP10')
       end if
    end if
!if(debug_flag) write(*,"(2a,4i5,2g12.3)") "GMO3 ", trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, spod, epod, L%FstO3,L%cano3_ppb

  end subroutine Calc_POD 

 !=========================================================================
  subroutine Calc_SPOD(iO3cl,iLC, spod, debug_flag, debug_txt )
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: spod

    real    :: o3, Y
    integer :: i,j
    real :: spod_fenv, fgmax
    real :: dg, bt, dTs, tmp, f_temp, f_vpd !OWN CALCS

    spod = 0.0

    if ( Grid%izen >= AOT_HORIZON ) RETURN  !UN or MM use daylight

   ! Start and end of POD accumulation period:

   ! Outside growing season. We check agsinst f_phen also, since for jday==SGS
   ! deciduous forests get f_phen zero
    if ( jday < L%SGS .or. jday > L%EGS .or. L%f_phen < 1.0e-5 ) then
       return
    end if

    i = Grid%i
    j = Grid%j
    Y = VEGO3_OUTPUTS(iO3cl)%Threshold ! nmole/m2/s

    o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) * cfac(IXADV_O3,i,j) * PPBINV

    fgmax = 1.0
    if ( index( VEGO3_OUTPUTS(iO3cl)%name, "spruce" )>0 ) fgmax = 0.57
   ! For crops, we use the simple Table 3.14 from MM, to get from 3m to 1m
    if ( index( VEGO3_OUTPUTS(iO3cl)%name, "crops" )>0 )  fgmax = 2.5 * 0.88/0.95

!..3) Calculate  f_temp
!---------------------------------------
! Asymmetric  function from Mapping Manual
! NB _ much more efficient to tabulate this - do later!
  
  dg  =    ( do3se(iLC)%T_opt - do3se(iLC)%T_min )
  bt  =    ( do3se(iLC)%T_max - do3se(iLC)%T_opt ) / dg
  dTs = max( do3se(iLC)%T_max - Grid%t2C, 0.0 )      !CHECK why max?
  tmp = dTs / ( do3se(iLC)%T_max - do3se(iLC)%T_opt )
  f_temp = ( Grid%t2C - do3se(iLC)%T_min ) / dg *  tmp**bt

  f_temp = max( L%f_temp, 0.01 )  ! Revised usage of min value during 2007


!..4) Calculate f_vpd
!---------------------------------------

 f_vpd = do3se(iLC)%f_min + &
          (1.0-do3se(iLC)%f_min) * (do3se(iLC)%VPD_min - L%vpd )/ &
              (do3se(iLC)%VPD_min - do3se(iLC)%VPD_max )
 f_vpd = min(f_vpd, 1.0)
 f_vpd = max(f_vpd, do3se(iLC)%f_min)


    if (DEBUG%AOT .and.  debug_flag ) then
      write( *,*) "spod_fenv POS"//trim(VEGO3_OUTPUTS(iO3cl)%name), me, i_fdom(i), j_fdom(j)
      write(*,"(a,i3,i4,2i3,f7.3,i5,2es10.2,10f8.3)") "spod_fenv:", iLC,  &
         jday, current_date%month, current_date%day, &
         current_date%hour + current_date%seconds/3600.0,& 
         Grid%izen, L%PARsun, Grid%Idirect,  &
         L%f_light, L%f_temp, L%rh,&
         L%f_vpd, L%f_phen
    end if
    spod_fenv = L%f_temp *  L%f_vpd  *  L%f_phen * fgmax
    spod_fenv = max( spod_fenv ,  L%f_min)

    spod = max ( o3 * spod_fenv - Y, 0.0 )

    if (DEBUG%AOT .and.  debug_flag ) then

      write(*,"(a,i4,2i3,f6.2,2i4,1x,3f7.2,i3,6f6.2,2x,2f7.2)") "OPOD"//trim(VEGO3_OUTPUTS(iO3cl)%name), &
         jday, current_date%month, current_date%day, &
         current_date%hour + current_date%seconds/3600.0,& 
         L%SGS, L%EGS,&
         o3, L%cano3_ppb, o3 * spod_fenv, nint(Y), &
              L%f_env,  spod_fenv, L%f_temp, L%f_vpd, &
               L%f_phen, L%f_light, spod, L%FstO3

      if ( present( debug_txt )) then
        write(6,"(a,L2,3i4)") "SPODACC ", VEGO3_OUTPUTS(iO3cl)%RelSGS, &
           VEGO3_OUTPUTS(iO3cl)%SAccPeriod, VEGO3_OUTPUTS(iO3cl)%EAccPeriod, &
           jday
        call datewrite("YYY"//trim(debug_txt) // "defn:" // &
          trim(VEGO3_OUTPUTS(iO3cl)%defn), iLC, &
           (/ real(Grid%izen), Y,  o3, spod /) )
      end if
    end if

  end subroutine Calc_SPOD 

 !=========================================================================


! Will move SOMO and maxo3 here in future versions
!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September


end module AOTx_ml
