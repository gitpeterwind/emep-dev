!_____________________________________________________________________________!
 module Setup_1d_ml

  ! DESCRIPTION
  ! Generates arrays for 1-D column , for inpout to chemical solver. The output
  ! fields are stored in the Setup_1dfields_ml module.


  !-----------------------------------------------------------------------!
  use Volcanos_ml
  use AirEmis_ml,            only :  airn, airlig   ! airborne NOx emissions
  use Biogenics_ml         , only :  emnat,canopy_ecf, BIO_ISOP, BIO_TERP
  use Chemfields_ml,         only :  xn_adv,xn_bgn,xn_shl         
  use Dates_ml,              only : date, dayno
  use Emissions_ml,          only :  gridrcemis, KEMISTOP    !hf VOL
  use Functions_ml,          only :  Tpot_2_T                !ds apr2005
  use GenSpec_tot_ml,        only :  SO4,aNO3,pNO3
  use GenSpec_adv_ml,        only :  NSPEC_ADV
  use GenSpec_shl_ml,        only :  NSPEC_SHL
  use GenSpec_bgn_ml,        only :  NSPEC_COL, NSPEC_BGN &
    ,xn_2d_bgn         !u2 concentration terms, specified species
  use MyChem_ml,             only :  Set_2dBgnd
  use GenRates_rct_ml,       only :  NRCT, rcit !u3 Tabulated rate coeffs
  use GenRates_rcmisc_ml,    only :  NRCMISC, set_rcmisc_rates
  use GridValues_ml,         only :  sigma_mid, xmd, carea, i_glob, j_glob
  use MassBudget_ml,         only :  totem    ! sum of emissions
  use Met_ml,                only :  roa, th, ps, q, t2_nwp, cc3dmax &
                                    ,zen, Idirect, Idiffuse &
                                    ,z_bnd     !ds Pb210
  use ModelConstants_ml,     only :  &
     ATWAIR                          &        
    ,dt_advec                        & ! time-step
    ,current_date                    & ! 
    ,PT                              & ! Pressure at top
    ,MFAC                            & ! converts roa (kg/m3 to M, molec/cm3)
    ,KMAX_MID ,KMAX_BND, KCHEMTOP                ! Start and upper k for 1d fields
  use My_Aerosols_ml,    only : SEASALT        !SeaS
  use My_Emis_ml,           only : NRCEMIS  , AIRNOX, QRCAIRNO &
                                  ,QRCAIRNO2, NFORESTVOC&
                                  ,QRCVOL,VOLCANOES &  ! hf -extended VOL
                                  ,NSS  !SeaS
  use My_MassBudget_ml,      only : N_MASS_EQVS, ixadv_eqv, qrc_eqv
  use My_BoundConditions_ml, only : BGN_2D     !hf u2
  use Par_ml,                only :  me& !!(me for tests)
                             ,gi0,gi1,gj0,gj1,ISMBEG,JSMBEG !hf VOL
  use PhysicalConstants_ml,  only :  AVOG, PI
  use Radiation_ml,          only : & !ds mar2005 zen, Idirectt, Idiffuse, 
                              PARfrac, Wm2_uE  ! ds rv1_6_x for bio
  use Setup_1dfields_ml,     only : &
     xn_2d                &  ! concentration terms
    ,rcemis, rcbio        &  ! emission terms
    ,rc_Rn222             &  ! ds Pb210
    ,rct, rcmisc          &  ! emission terms
    ,rcss                 &  !SeaS - sea salt
    ,rh, temp, itemp,pp      &  ! 
    ,amk                  &  ! Air concentrations 
    ,izen &
    ,f_Riemer  !weighting factor for N2O5 hydrolysis    
                     ! integer of zenith angle
   use SeaSalt_ml,        only : SS_prod   !SeaS
   use Tabulations_ml,    only :  tab_esat_Pa
   use UKdep_ml,          only : water_fraction, ice_fraction  !dsPb210
  implicit none
  private
  !-----------------------------------------------------------------------!


  public :: setup_1d   ! Extracts results for i,j column from 3-D fields
  public :: setup_bio    ! Biogenic emissions
  public :: setup_rcemis ! Emissions  (formerly "poll")
  public :: reset_3d     ! Exports results for i,j column to 3-D fields


contains
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   subroutine setup_1d(i,j)

 !..   extracts data along one vertical column for input to chemical
 !     solver concentrations for chemistry......
 !
 ! Outputs, amk, o2k, rcairlig, ...

    integer, intent(in) :: i,j    ! coordinates of column

   !/* local

    integer           :: k, n, ispec, irc    ! loop variables
    real              :: qsat ! saturation water content


    real ,dimension(KCHEMTOP:KMAX_MID) :: tinv, & ! Inverse of temp.
        h2o, o2k     ! water, O2

    ! - calculate integer of local zenith angle here

     izen = int ( zen(i,j) + 0.5 )
     izen = max(1, izen)    ! Just to avoid zero in indices.

    do k = KCHEMTOP, KMAX_MID
 
  !- MFAC replaces earlier use of CHEFAC and ATWAIR - to scale from
  !  density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

       amk(k) = roa(i,j,k,1) * MFAC  ! molecules air/cm3


       h2o(k) = max( 1.e-5*amk(k), &
                     q(i,j,k,1)*amk(k)*ATWAIR/18.0)

      ! nb. max function for h2o  used as semi-lagrangian scheme used
      ! in LAM50 (and HIRLAM) often gives negative H2O....   :-(


       pp(k) = PT + sigma_mid(k)*(ps(i,j,1) - PT)

       !ds apr2005 temp(k) = th(i,j,k,1)*exp(XKAP*log(pp(k)*1.e-5))
       temp(k) = th(i,j,k,1)* Tpot_2_T( pp(k) )

       itemp(k) = nint( temp(k) -1.E-9)
!pw the "-1.E-9" is put in order to avoid possible different roundings on different machines. 

       ! relative humidity - moved here from setup_miscrc. Note that the
       ! MADE/MACHO codes used sigma_bnd instead of sigma_mid for this pp, but
       ! sigma_mid seems better.

       !old pp   = PT + sigma_bnd(k)*(ps(i,j,1) - PT)

       qsat  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
       rh(k) = min( q(i,j,k,1)/qsat , 1.0) 

        ! 1)/ Short-lived species - no need to scale with M
        !  13/9/2002 - reintroduced reset of xn_2d

         do n = 1, NSPEC_SHL
               xn_2d(n,k) = max(0.0,xn_shl(n,i,j,k))
         end do ! ispec

        ! 2)/ Advected species
        do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_2d(ispec,k) = max(0.0,xn_adv(n,i,j,k)*amk(k))
        end do ! ispec
 
        ! 3)/ Background species
        !** CTM2 with units in mix. ratio
        do n = 1, NSPEC_BGN
              xn_2d_bgn(n,k) = max(0.0,xn_bgn(n,i,j,k)*amk(k))
        end do ! ispec   

!hf setup weighting factor for hydrolysis  
             f_Riemer(k)=96.*xn_2d(SO4,k)/( (96.*xn_2d(SO4,k))+(62.*xn_2d(aNO3,k)) )

   end do ! k


   o2k(:) = 0.21*amk(:)
   tinv(:) = 1./temp(:)

  !u2 - NEW ****
  !u3 - amk added to arguments

!hf   if (BGN_2D) &
!hf     call Set_2dBgnd(izen,cc3dmax(i,j,:),amk)  ! Sets remaining xn_2d_bgn

  ! 5 ) Rates 

   rct(:,:)    = rcit(:,itemp(:))

   !u1 call set_rctroe_rates(tinv,amk,rctroe)

   call set_rcmisc_rates(itemp,tinv,amk,o2k,h2o,rh,rcmisc)




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
     integer ::  iqrc,k,n
     real    :: scaling, scaling_k
     real    :: eland   !ds Pb210  - emissions from land

!hf VOL
    integer ::  i_help,j_help,i_l,j_l

!hf VOL initilize
    rcemis(:,:)=0.    
    rcss(:,:) = 0.  !SeaS 
   
     !rv1.2.1 do k=KMAX_MID-3,KMAX_MID -- bug spotted by st, 2/12/2002
     do k=KEMISTOP,KMAX_MID

        do iqrc = 1, NRCEMIS
          rcemis(iqrc,k) = gridrcemis(iqrc,k,i,j)
        end do ! iqrc
     enddo

     !/** Add volcanoe emissions
    
     if ( VOLCANOES  ) then ! for models that include volcanos
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
              rcemis(QRCVOL,k)=rcemis(QRCVOL,k)+rcemis_volc(volc_no)
              !write(*,*)'Adding volc. emissions ',rcemis_volc(volc_no),volc_no,&
              !            'to height=',k,'i,j',i_help,j_help
              !write(*,*)'TOT rcemis=',rcemis(QRCVOL,:)
           endif
        endif
     enddo
    
     endif ! VOLCANOES

    !/** lightning and aircraft ... Airial NOx emissions if required:

     if ( AIRNOX  ) then

       !QRCAIRNO is set to QRCNO and QRCAIRNO2 is set to QRCNO2 if 
       ! AIRNOX is true. Otherwise to a
       ! dummy value of 1. Avoids problems with
       !undefined QRCNO in non-NOx models.

        !rv1.2.1 do k=KCHEMTOP,KMAX_MID-4  ! bug, 2/12/2002
        do k=KCHEMTOP,KEMISTOP-1
          rcemis(QRCAIRNO,k)  = 0.95 * (airn(k,i,j)+airlig(k,i,j))
          rcemis(QRCAIRNO2,k) = 0.05 * (airn(k,i,j)+airlig(k,i,j))

        enddo

        !rv1.2.1 do k=KMAX_MID-3,KMAX_MID  ! bug, 2/12/2002
        do k=KEMISTOP,KMAX_MID
          rcemis(QRCAIRNO,k)  = rcemis(QRCAIRNO,k) &  
                              + 0.95 * (airn(k,i,j)+airlig(k,i,j))
          rcemis(QRCAIRNO2,k) = rcemis(QRCAIRNO2,k) &
                              + 0.05 * (airn(k,i,j)+airlig(k,i,j))

        enddo

     end if ! AIRNOX

!SeaS......
     !/** Add sea salt production
    
     if ( SEASALT  ) then

          do iqrc = 1, NSS
            rcss(iqrc,KMAX_MID) = SS_prod(iqrc,i,j)
          enddo

     endif
!SeaS......

!6b Mass Budget calculations
!  jej Adding up the emissions in each timestep
!  ds  use ixadv_eqv, qrc_eqv from My_Emis_ml, e.g. ixadv_eqv(1) = IXADV_SO2,
!      qrc_eqv(1) = QRCSO2

   scaling = dt_advec * xmd(i,j) * (ps(i,j,1) - PT)

   do k = KCHEMTOP,KMAX_MID

       scaling_k = scaling * carea(k)/amk(k)

       do n = 1, N_MASS_EQVS  
          totem( ixadv_eqv(n) ) = totem( ixadv_eqv(n) ) + &
                rcemis( qrc_eqv(n),k) * scaling_k
       end do

   end do ! k loop 

  ! Soil Rn222 emissions from non-ice covered land, + water
  ! at rate of 1 atom/cm2/s

     eland = 1.0 - water_fraction(i,j) - ice_fraction(i,j)
!hf initialize, needed in My_Reactions
     rc_Rn222(:)=0.0     
!     rc_Rn222(KMAX_MID) = &
!            ( 0.00182 * water_fraction(i,j)  + eland ) / &
!            (z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND)) 
!hf z_bnd is in m, not cm, so need to divide by 100.
     rc_Rn222(KMAX_MID) = &
            ( 0.00182 * water_fraction(i,j)  + eland ) / &
            ((z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND))*100.) 

  end subroutine setup_rcemis
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine setup_bio(i,j)
  !
  !---- assign isoprene rates  ------------------------------------------------
  !
  !  So far, assigns isoprene using surface (2m) temperature, and for all
  !  zenith angles <90. Should include light dependance at some stage
  !
  !  Output : rcbio - isoprene emissions for 1d column
  !
  !  Called from setup_ml, every  advection step.
  !----------------------------------------------------------------------------

!  input
  integer, intent(in) ::  i,j

!  local
  integer la,it2m,n,k,base,top,iclcat
  real clear

 !ds rv1_6_x  - light effects added for isoprene emissions

  real            :: par   !ds Photosynthetically active radiation
  real            :: cL    !ds Factor for light effects
  real, parameter :: &
      CL1 = 1.066  , &    !ds Alex Guenther's params
      ALPHA = 0.0027      !ds Alex Guenther's params

  if ( NFORESTVOC == 0  ) return   ! e.g. for MADE



  it2m = nint(t2_nwp(i,j,1)-273.15-1.E-9)
!pw the "-1.E-9" is put in order to avoid possible different roundings on different machines.
  it2m = max(it2m,1)
  it2m = min(it2m,40)

  rcbio(BIO_TERP,KMAX_MID) = emnat(BIO_TERP,i,j)*canopy_ecf(BIO_TERP,it2m)

  ! Isoprene has emissions in daytime only:
  rcbio(BIO_ISOP,:) = 0.0
  if (izen <= 90) then

     ! ds light effects from Guenther G93
      par = (Idirect(i,j) + Idiffuse(i,j)) * PARfrac * Wm2_uE
      cL = ALPHA * CL1 * par/ sqrt( 1 + ALPHA*ALPHA * par*par)

      rcbio(BIO_ISOP,KMAX_MID) = emnat(BIO_ISOP,i,j)*canopy_ecf(BIO_ISOP,it2m) * cL
  endif

  end subroutine setup_bio

  !----------------------------------------------------------------------------
   subroutine reset_3d(i,j)

    integer, intent(in) :: i,j
    integer :: k, n, ispec, irc    ! loop variables
 
 
         do k = KCHEMTOP, KMAX_MID
 

           ! 1)/ Short-lived species - no need to scale with M
           !  13/9/2002 - reintroduced reset of xn_2d

            do n = 1, NSPEC_SHL
               xn_shl(n,i,j,k) = xn_2d(n,k)
            end do ! ispec
 
           ! 2)/ Advected species

           do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
!              xn_adv(n,i,j,k) = max(0.0, xn_2d(ispec,k)/amk(k))
              xn_adv(n,i,j,k) = xn_2d(ispec,k)/amk(k)
           end do ! ispec
  
         end do ! k
 
   end subroutine reset_3d
   !---------------------------------------------------------------------------

end module Setup_1d_ml
!_____________________________________________________________________________!













