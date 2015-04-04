module Setup_1d_ml
!-----------------------------------------------------------------------!
! DESCRIPTION
! Generates arrays for 1-D column , for input to chemical solver. The output
! fields are stored in the Setup_1dfields_ml module.
!-----------------------------------------------------------------------!

use AeroFunctions,       only: umWetRad, pmSurfArea, cMolSpeed, UptakeRate
use AirEmis_ml,          only: airn, airlig   ! airborne NOx emissions
use Biogenics_ml,        only: SoilNOx
use Biogenics_ml,        only: EMIS_BioNat, EmisNat  
use Chemfields_ml,       only: xn_adv,xn_bgn,xn_shl, &
                               NSPEC_COL, NSPEC_BGN, xn_2d_bgn
use ChemGroups_ml,       only: WDEP_PM10_GROUP, chemgroups
use CheckStop_ml,        only:  CheckStop
use DerivedFields_ml,    only: d_2d
use EmisGet_ml,          only:  nrcemis, iqrc2itot  !DSRC added nrcemis
use Emissions_ml,        only:  gridrcemis, gridrcroadd, SumSplitEmis, KEMISTOP
use ForestFire_ml,       only: Fire_rcemis, burning
use Functions_ml,        only:  Tpot_2_T
use ChemFields_ml,       only: SurfArea_um2cm3
use ChemSpecs  !,           only:  SO4,C5H8,NO,NO2,SO2,CO,
use ChemRates_rct_ml,    only:  set_rct_rates, rct
use GridValues_ml,       only:  xmd, GridArea_m2, & 
                                 debug_proc, debug_li, debug_lj,&
                                 A_mid,B_mid,gridwidth_m,dA,dB,&
                                 i_fdom, j_fdom
use Io_Progs_ml,         only: datewrite !MASS
use LocalVariables_ml,   only: Grid
use MassBudget_ml,       only: totem    ! sum of emissions
use MetFields_ml,        only: ps
use MetFields_ml,        only: roa, th, q, t2_nwp, cc3dmax, &
                               zen, Idirect, Idiffuse,z_bnd
use ModelConstants_ml,   only:  &
   DEBUG,DEBUG_MASS             &
!EXCL    ,USES, MINCONC                & ! conc.limit if USES%MINCONC
  ,AERO                         & ! for wet radii and surf area.
  ,dt_advec                     & ! time-step
  ,PT                           & ! Pressure at top
  ,USES                         & ! Forest fires so far
  ,USE_SEASALT                  &
  ,USE_LIGHTNING_EMIS, USE_AIRCRAFT_EMIS      & !
  ,USE_GLOBAL_SOILNOX, USE_DUST, USE_ROADDUST & !
  ,USE_EMERGENCY,DEBUG_EMERGENCY   & ! Emergency: Volcanic Eruption
  ,KMAX_MID ,KMAX_BND, KCHEMTOP      ! Start and upper k for 1d fields
use My_Derived_ml,       only: EmisSplit_OUT
use Landuse_ml,          only: water_fraction, ice_landcover
use Par_ml,              only: me,MAXLIMAX,MAXLJMAX, & 
                               gi0,gi1,gj0,gj1,IRUNBEG,JRUNBEG
use PhysicalConstants_ml,only: ATWAIR, AVOG, PI, GRAV
use Radiation_ml,        only: PARfrac, Wm2_uE
use Setup_1dfields_ml,   only: &
   xn_2d                &  ! concentration terms
  ,rcemis               &  ! emission terms
  ,rh, temp, tinv, itemp,pp      &  !
  ,amk, o2, n2, h2o     &  ! Air concentrations
  ,cN2O5, cHO2, cO3, cHNO3 &  ! mol speeds, m/s
  ,ugdryPM,DpgNw,S_m2m3    ! for wet diameter and surf area
use SmallUtils_ml,       only: find_index
use Tabulations_ml,      only: tab_esat_Pa
use TimeDate_ml,         only: current_date, date
use Volcanos_ml
use Emergency_ml,        only: EmergencyRate
use Units_ml,            only: to_number_cm3 ! converts roa [kg/m3] to M [molec/cm3]

implicit none
private
!-----------------------------------------------------------------------!

public :: setup_1d   ! Extracts results for i,j column from 3-D fields
public :: setup_rcemis ! Emissions  (formerly "poll")
public :: reset_3d     ! Exports results for i,j column to 3-D fields

! Indices for the species defined in this routine. Only set if found
! Hard-coded for 2 specs just now. Could extend and allocate.
integer, private, parameter :: NROADDUST = 2
integer, private, parameter :: iROADF=1,  iROADC=2
integer, private, save :: inat_RDF,  inat_RDC, inat_Rn222
integer, private, save :: itot_RDF=-999,  itot_RDC=-999, itot_Rn222=-999

! Minimum concentration allowed for advected species. Avoids some random
! variations between debugged/normal runs with concs of e.g. 1.0e-23 ppb or
! less. 1.0e-20 ppb is ca. 2.5e-10 molec/cm3 at sea-level, ca. 2.5e-11 at top
! e.g.: MINCONC = 1.0e-29

!DUST_ROAD_F

contains
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   subroutine setup_1d(i,j)

 !..   extracts data along one vertical column for input to chemical
 !     solver concentrations for chemistry......
 !

    integer, intent(in) :: i,j    ! coordinates of column
    character(len=9)  :: sub='setup_1d:'
    character(len=30)  :: fmt="(a,i3,99g12.3)"  ! default format
    logical :: debug_flag
    logical, save :: first_call = .true.
    real :: wradR,wradS,wradU, ugdryPM, ugDustF, ugSSaltF
    real, dimension(6)  :: umRdry ! 6 is hard-coded equiv of AERO%NSAREA.
    real :: wradSC,wradUC, ugDustC, ugSSaltC
    integer :: iw, ipm ! for wet rad

   !** local

    integer           :: k, n, ispec    ! loop variables
    real              :: qsat ! saturation water content

    debug_flag =  ( DEBUG%SETUP_1DCHEM .and. debug_proc .and.  &
      i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )
    if( debug_flag ) write(*,*) sub//"====================  "


    do k = KCHEMTOP, KMAX_MID

  !- to_number_cm3 - to scale from  density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

       amk(k) = roa(i,j,k,1) * to_number_cm3  ! molecules air/cm3

       h2o(k) = max( 1.e-5*amk(k), &
                     q(i,j,k,1)*amk(k)*ATWAIR/18.0)

      ! nb. max function for h2o  used as semi-lagrangian scheme used
      ! in LAM50 (and HIRLAM) often gives negative H2O....   :-(

       pp(k) = A_mid(k) + B_mid(k)*ps(i,j,1)

       temp(k) = th(i,j,k,1)* Tpot_2_T( pp(k) )

       itemp(k) = nint( temp(k) -1.E-9) ! the "-1.E-9" is put in order to
               ! avoid possible different roundings on different machines.

       qsat  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
       rh(k) = min( q(i,j,k,1)/qsat , 1.0)
       rh(k) = max( rh(k) , 0.001)

        ! 1)/ Short-lived species - no need to scale with M

         do n = 1, NSPEC_SHL
               xn_2d(n,k) = max(0.0,xn_shl(n,i,j,k))
         end do ! ispec

        ! 2)/ Advected species
        do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_2d(ispec,k) = max(0.0,xn_adv(n,i,j,k)*amk(k))
!EXp.             ! Avoids tiny numbers that complicate comparisons in DEBUG mode
!EXp.              if( USES%MINCONC ) then
!EXp.                if( xn_2d(ispec,k) < MINCONC ) xn_2d(ispec,k) = 0.0
!EXp.              end if
        end do ! ispec

        ! 3)/ Background species ( * CTM2 with units in mix. ratio)
        do n = 1, NSPEC_BGN
              xn_2d_bgn(n,k) = max(0.0,xn_bgn(n,i,j,k)*amk(k))
        end do ! ispec

      ! Surf Area
        if ( USES%SURF_AREA ) then ! GERBER

           ispec=NO3_c ! CRUDE HARD CODE for now
           ugdryPM     = 0.27*xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG

           ugDustF  = 0.0
           ugSSaltF = 0.0
           ugDustC  = 0.0
           ugSSaltC = 0.0

           do ipm = 1, size( WDEP_PM10_GROUP ) !! PMFINE_GROUP )
             ispec = WDEP_PM10_GROUP(ipm)

             if( index( species(ispec)%name, 'SEASALT_F' )>0) then

                ugSSaltF = ugSSaltF + xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
           if( debug_flag .and. k==20 )  print fmt, "GERB ugSS "// trim(species(ispec)%name),ispec, ugSSaltF
             else if( index( species(ispec)%name, 'SEASALT_C' )>0) then
                ugSSaltC = ugSSaltC + xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG

             else if ( index( species(ispec)%name, 'DUST' )>0) then
                if( index( species(ispec)%name, '_F' )>0) then
                    ugDustF = ugDustF + xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG

           if( debug_flag .and. k==20 )  print fmt, "GERB ugDU "// trim(species(ispec)%name),ispec, ugDUSTF
                else
                    ugDustC = ugDustC + xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
                end if

             else ! don't distinguish F from C for NSD
               ugdryPM = ugdryPM + &
                xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG !-> ug/m3
           if( debug_flag .and. k==20 )  print fmt, "GERB ugPM "// trim(species(ispec)%name),ispec, ugdryPM
             end if
           end do !ipm

          ! call, for fine aeroso
           ! excludign dust and sea-salt,
           ! GERBER equations for wet radius

          DpgNw(AERO%NSD_F,:)=AERO%DpgN(1)
          DpgNw(AERO%NSD_C,:)=AERO%DpgN(2)
          DpgNw(AERO%SS_F,:)=AERO%DpgN(1)
          DpgNw(AERO%SS_C,:)=AERO%DpgN(3)
          DpgNw(AERO%DU_F,:)=AERO%DpgN(1)
          DpgNw(AERO%DU_C,:)=AERO%DpgN(4)

          do iw = 1, AERO%NSAREA 
           umRdry(iw) =  0.5e6*AERO%DpgN( AERO%Ddry(iw))
           if( AERO%Gb(iw) > 0 ) then
              DpgNw(iw,k)  = 2 * 1.0e-6*umWetRad( umRdry(iw), rh(k), AERO%Gb(iw)) 
           else ! index -1 indicates use dry (for dust)
              DpgNw(iw,k)  = AERO%DpgN( AERO%Ddry(iw))
           end if
          end do

           iw= AERO%NSD_F
           S_m2m3(iw,k) = pmSurfArea(ugdryPM,Dp=2*umRdry(iw), Dpw=DpgNw(iw,k))

           iw= AERO%SS_F
           S_m2m3(iw,k) = pmSurfArea(ugSsaltF,Dp=2*umRdry(iw), Dpw=DpgNw(iw,k))

           iw= AERO%SS_C
           S_m2m3(iw,k) = pmSurfArea(ugSsaltC,Dp=2*umRdry(iw), Dpw=DpgNw(iw,k))

          ! dust - just use dry radius
           iw= AERO%DU_F
           S_m2m3(iw,k) = pmSurfArea(ugDustF,Dp=2*umRdry(iw), Dpw=2*umRdry(iw))
           iw= AERO%DU_C
           S_m2m3(iw,k) = pmSurfArea(ugDustC,Dp=2*umRdry(iw), Dpw=2*umRdry(iw))

!           call CheckStop (  S_m2m3(k) < 0.0 , "NEGS_m2m3" )
!           if( debug_flag .and. k==20 )  print fmt, "GERB ugdry ", k, wradR*1.0e6, S_m2m3(k)*1.0e6

           if( debug_flag .and. k==20 )  then
            print fmt,  sub//" SAREAugPM in  ", k, ugdryPM, ugSSaltC, ugDustC
            do iw = 1, AERO%NSAREA
             print fmt, sub//"GERB ugDU  ", iw, umRdry(iw), S_m2m3(iw,k)
            end do
           end if

           S_m2m3(:,k) = min( S_m2m3(:,k), 6.0e-3)  !! Allow max 6000 um2/cm3

           ! m2/m3 -> um2/cm3 = 1.0e6, for output to netcdf
           if( k == KMAX_MID ) then 
              do iw = 1, AERO%NSAREA
                SurfArea_um2cm3(iw,i,j) = 1.0e6* S_m2m3(iw,k)
              end do
           end if

!         if( me==36 .and. i==5.and.j==31.and. Sarea_um2(kpm) > 1.0e5 ) then
!             print "(a,4i4,g12.3,f7.2,es12.3)", "SURFAREA ",me,i,j,kpm,
!             ugdryPM(kpm), rh(kpm), Sarea_um2(kpm)
!             !call StopAll('SSS')
!         end if

         !! Mol. speeds for aerosol uptake

          cn2o5 = cMolSpeed(temp,108.0)
          chno3 = cMolSpeed(temp, 63.0)
          cho2  = cMolSpeed(temp, 33.0)
          co3   = cMolSpeed(temp, 48.0)

! if( debug_flag .and. k==20 ) then
!    print fmt, "SSGAMMA ",k, rh(k), rct(71,k),cN2O5(k), DpgNWS(k),S_SSFM2M3(k),
!    UPTAKERATE(0.5*DPGNWS(k),CN2O5(k),gam=GammaN2O5_EJSS(rh(k)),S=S_SSFM2M3(k))
!end if

         end if ! GERBER


      ! End Surf Area
    
   end do ! k

! Check that concentrations are not "contaminated" with NaN
   if ( isnan( xn_2d(NO2,KMAX_MID) ) ) then
      print *, "NANAN ", trim(species(NO2)%name)
      call CheckStop( "Detected non numerical concentrations (NaN)")
   end if

   o2(:) = 0.21 *amk(:)
   n2(:) = amk(:) - o2(:)
!   o2(:) = 0.2095 *amk(:) ! more exact, but prefer o3+n2 to add to 100%
!   n2(:) = 0.7808 *amk(:)
   tinv(:) = 1./temp(:)



  ! 5 ) Rates  (!!!!!!!!!! NEEDS TO BE AFTER RH, XN, etc. !!!!!!!!!!)


   call set_rct_rates()

   if ( first_call ) then
     call CheckStop( any(isnan(rct(:,:))), sub//"RCT NAN'd")
     first_call = .false.
   end if

  if( debug_flag ) then
      write(*,"(a,10es10.3)") sub//"RCT ", &
            rct(3,KMAX_MID), rct(4,KMAX_MID)
      write(*,"(a,10es10.3)") sub//"XN  ", &
        amk(KMAX_MID),  xn_2d(O3,KMAX_MID),  xn_2d(NO2,KMAX_MID)
      write(*,"(a,10es10.3)") sub//"1D-Riemer", xn_2d(SO4,KMAX_MID),&
        xn_2d(NO3_F,KMAX_MID)
  end if



   end subroutine setup_1d
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine setup_rcemis(i,j)

    !-------------------------------------------------------------------
    !**    DESCRIPTION:
    !    Extracts emissions in column from gridrcemis, for input to chemistry
    !    routines. Results in "rcemis" array
    !-------------------------------------------------------------------

   !-- arguments
     integer, intent(in) ::  i,j     ! coordinates of column

   !  local
     integer ::  iqrc,k, itot
     real    :: eland   ! for Pb210  - emissions from land

    integer ::  i_help,j_help,i_l,j_l
    logical, save     :: first_call = .true. 

    if ( first_call ) then
      inat_RDF = find_index( "DUST_ROAD_F", EMIS_BioNat(:) )
      inat_RDC = find_index( "DUST_ROAD_C", EMIS_BioNat(:) )
      itot_RDF = find_index( "DUST_ROAD_F", species(:)%name    )
      itot_RDC = find_index( "DUST_ROAD_C", species(:)%name    )
      itot_Rn222= find_index( "RN222", species(:)%name    )
      first_call = .false.
    end if 

! initilize ! initilize ! initilize ! initilize
    rcemis(:,:)=0.
! initilize ! initilize ! initilize ! initilize

     do k=KEMISTOP,KMAX_MID

        do iqrc = 1, NRCEMIS
          itot = iqrc2itot( iqrc )
          rcemis(itot,k) = gridrcemis( iqrc ,k,i,j)
        end do ! iqrc
     enddo

     !*** Add volcanoe emissions

     if ( Volcanoes_found  ) then ! for models that include volcanos
                      !QRCVOL=QRCSO2 for models with volcanos
                      !For non-volc models it is set to dummy value 1
                      !to avoid problems with undefined QRCVOL
     do volc_no=1,nvolc
        i_help=i_volc(volc_no) !Global coordinates
        j_help=j_volc(volc_no) !Global coordinates
        if ((i_help >= gi0).and.(i_help<=gi1).and.(j_help>= gj0)&
             .and.(j_help<= gj1))then !on the correct processor
           i_l=i_help - gi0  +1
           j_l=j_help - gj0  +1
           if((i_l==i).and.(j_l==j))then !i,j have a volcano
              k=height_volc(volc_no)
              rcemis(SO2,k)=rcemis(SO2,k)+rcemis_volc(volc_no)
              !write(*,*)'Adding volc. emissions ' &
              !          ,rcemis_volc(volc_no),volc_no,&
              !          'to height=',k,'i,j',i_help,j_help
              !write(*,*)'TOT rcemis=',rcemis(QRCVOL,:)
           endif
        endif
     enddo

     endif ! VOLCANOES

! Contribution from Emergeny scenarios:
!   Volcanic eruption scenarios set SO2 and ASH tracers emission rates.
!   NPP accident scenarios set redoactive tracers emission rates.
!   EmergencyRate(i,j) returns an array of rcemis size.
    if(USE_EMERGENCY) rcemis(:,:)=rcemis(:,:)+EmergencyRate(i,j)

    !*** lightning and aircraft ... Airial NOx emissions if required:

     if ( USE_LIGHTNING_EMIS  ) then

        do k=KCHEMTOP, KMAX_MID

          rcemis(NO,k)  = rcemis(NO,k) &
                              + 0.95 * airlig(k,i,j)
          rcemis(NO2,k) = rcemis(NO2,k) &
                              + 0.05 * airlig(k,i,j)

        enddo
        if ( DEBUG%SETUP_1DCHEM .and. debug_proc .and.  &
               i==debug_li .and. j==debug_lj ) write(*,"(a,10es10.3)") &
                 " DEBUG_SETUP_AIRNOX ", airn(KMAX_MID,i,j),airlig(KMAX_MID,i,j)

     end if ! LIGHTNING_EMIS
     if ( USE_AIRCRAFT_EMIS  ) then

        do k=KCHEMTOP, KMAX_MID

          rcemis(NO,k)  = rcemis(NO,k) &
                              + 0.95 * airn(k,i,j)
          rcemis(NO2,k) = rcemis(NO2,k) &
                              + 0.05 * airn(k,i,j)

        enddo
        if ( DEBUG%SETUP_1DCHEM .and. debug_proc .and.  &
               i==debug_li .and. j==debug_lj ) write(*,"(a,10es10.3)") &
                 " DEBUG_SETUP_AIRNOX ", airn(KMAX_MID,i,j),airlig(KMAX_MID,i,j)

     end if ! AIRCRAFT NOX

     !*** Add sea salt production

     if ( USE_ROADDUST .and. itot_RDF>0  ) then  ! Hard-code indices for now

        rcemis(itot_RDF,KMAX_MID) = gridrcroadd(1,i,j)
        rcemis(itot_RDC,KMAX_MID) = gridrcroadd(2,i,j)
     endif


     if ( USES%FOREST_FIRES   ) then
     if ( burning(i,j)  ) then

       call Fire_rcemis(i,j)

     endif 
     endif  !ForestFires

   !Soil NOx
     if( USE_GLOBAL_SOILNOX)then !NEEDS CHECKING NOV2011
        rcemis(NO,KMAX_MID)=rcemis(NO,KMAX_MID)+SoilNOx(i,j)
     endif

    if(EmisSplit_OUT)then
!put all added emissions in EmisSplit_OUT, also natural emissions
       SumSplitEmis(i,j,:)=0.0
      do k=KCHEMTOP, KMAX_MID
         do iqrc=1,nrcemis
            !give unit mg/m2/s dt_advec multiplied in Derived_ml 
            SumSplitEmis(i,j,iqrc) = SumSplitEmis(i,j,iqrc)+&
                 rcemis(iqrc2itot(iqrc),k)*species(iqrc2itot(iqrc))%molwt*&
                 (dA(k)+dB(k)*ps(i,j,1))/(GRAV*amk(k)*ATWAIR)
         enddo
      enddo
    endif


  ! Soil Rn222 emissions from non-ice covered land, + water
  ! at rate of 1 atom/cm2/s

     eland = 1.0 - water_fraction(i,j) - ice_landcover(i,j)

! z_bnd is in m, not cm, so need to divide by 100.

    if ( itot_Rn222 > 0 ) rcemis( itot_Rn222,KMAX_MID )  = & 
            ( 0.00182 * water_fraction(i,j)  + eland ) / &
            ((z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND))*100.)

!ESX     rc_Rnwater(KMAX_MID) = water_fraction(i,j)  / &
!ESX            ((z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND))*100.)

  end subroutine setup_rcemis

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine reset_3d(i,j)

    integer, intent(in) :: i,j
    integer :: k, n, ispec    ! loop variables


         do k = KCHEMTOP, KMAX_MID

           ! 1)/ Short-lived species - no need to scale with M

            do n = 1, NSPEC_SHL
               xn_shl(n,i,j,k) = xn_2d(n,k)
            end do ! ispec

           ! 2)/ Advected species

           do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_adv(n,i,j,k) = xn_2d(ispec,k)/amk(k)
           end do ! ispec

         end do ! k

   end subroutine reset_3d
   !---------------------------------------------------------------------------

end module Setup_1d_ml
!_____________________________________________________________________________!



!Experimental code. Will re-instate in future
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   !FUTURE  subroutine setup_nh3(i,j)  ! EXPERIMENTAL
   !FUTURE  !
   !FUTURE  !---- assign nh3 rates  ------------------------------------------------
   !FUTURE  !
   !FUTURE  !
   !FUTURE  !  use NH3 emission potential for each activity sector, T2 and wind
   !FUTURE  !
   !FUTURE  !  Output : rcnh3 - nh3 emissions for 1d column
   !FUTURE  !
   !FUTURE  !  Called from Setup_1d_ml, every advection timestep
   !FUTURE  !----------------------------------------------------------------------------
   !FUTURE    !  input
   !FUTURE    use calc_emis_potential_ml, only :emnh3
   !FUTURE    integer, intent(in) ::  i,j
   !FUTURE    real ::scaling,scaling_k
   !FUTURE    real :: snh3
   !FUTURE    integer :: k
   !FUTURE
   !FUTURE    !  local
   !FUTURE
   !FUTURE    rcnh3(:)=0.0
   !FUTURE    snh3=sum(emnh3(:,i,j))
   !FUTURE
   !FUTURE    if (NH3EMIS_VAR  )then
   !FUTURE
   !FUTURE       rcnh3(KMAX_MID) = snh3
   !FUTURE       scaling = dt_advec * xmd(i,j)* gridwidth_m*gridwidth_m / GRAV
   !FUTURE
   !FUTURE       do k = KCHEMTOP,KMAX_MID
   !FUTURE         scaling_k = scaling * (dA(k) + dB(k)*ps(i,j,1))/amk(k)
   !FUTURE
   !FUTURE          totem( IXADV_NH3 ) = totem( IXADV_NH3 ) + &
   !FUTURE               rcnh3(k) * scaling_k
   !FUTURE       enddo
   !FUTURE
   !FUTURE       if ( DEBUG_NH3 .and. i_fdom(i)==DEBUG_i .and. j_fdom(j)==DEBUG_j)then
   !FUTURE         write(6,*)'Tange coordinates, proc, i,j, ',me,i,j
   !FUTURE         write(6,*)'rcnh3',rcnh3(KMAX_MID)
   !FUTURE         write(6,*)'sum emnh3',sum(emnh3(:,i,j)),snh3
   !FUTURE         write(6,*)'emnh3 1',emnh3(1,i,j)
   !FUTURE       endif
   !FUTURE    else
   !FUTURE       rcnh3(KMAX_MID) =0.0
   !FUTURE    endif
