! <Setup_1d_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________!
 module Setup_1d_ml

  ! DESCRIPTION
  ! Generates arrays for 1-D column , for input to chemical solver. The output
  ! fields are stored in the Setup_1dfields_ml module.


  !-----------------------------------------------------------------------!
  use AirEmis_ml,            only :  airn, airlig   ! airborne NOx emissions
  use Biogenics_ml         , only :  emnat,canopy_ecf, BIO_ISOP, BIO_TERP
  use BoundaryConditions_ml, only : BGN_2D
  use Chemfields_ml,         only :  xn_adv,xn_bgn,xn_shl         
  use CheckStop_ml,          only :  CheckStop
  use EmisDef_ml,           only : AIRNOX, NBVOC,VOLCANOES, FOREST_FIRES &
                                  ,NSS  !SeaS
  use EmisGet_ml,          only :  nrcemis, iqrc2itot  !DSRC added nrcemis
  use Emissions_ml,          only :  gridrcemis, KEMISTOP
  use ForestFire_ml,         only : Fire_rcemis, burning !,  rcNOx_fire,rcCO_fire, rcHC_fire !DSFF
  use Functions_ml,          only :  Tpot_2_T
use ChemChemicals_ml,        only :  species
  use ChemSpecs_tot_ml,        only :  SO4,aNO3,pNO3,C5H8,NO,NO2,SO2,CO
  use ChemSpecs_adv_ml,        only :  NSPEC_ADV, IXADV_NO2, IXADV_O3
  use ChemSpecs_shl_ml,        only :  NSPEC_SHL
  use ChemSpecs_bgn_ml,        only :  NSPEC_COL, NSPEC_BGN, xn_2d_bgn
  use ChemRates_rct_ml,       only :  set_rct_rates, rct
  use ChemRates_rcmisc_ml,    only :  rcmisc, set_rcmisc_rates
  use GridValues_ml,         only :  sigma_mid, xmd, &
                                     debug_proc, debug_li, debug_lj,&
                                     A_mid,B_mid,gridwidth_m,dA,dB
  use LocalVariables_ml,     only :  Grid
  use MassBudget_ml,         only :  totem    ! sum of emissions
  use MetFields_ml,          only :  ps
  use MetFields_ml,                only :  roa, th, q, t2_nwp, cc3dmax &
                                    ,zen, Idirect, Idiffuse,z_bnd
  use ModelConstants_ml,     only :  &
     ATWAIR                          &        
    ,DEBUG_SETUP_1DCHEM              & !DSGC
    ,DEBUG_SETUP_1DBIO               & !DSGC
    ,dt_advec                        & ! time-step
    ,PT                              & ! Pressure at top
    ,MFAC                            & ! converts roa (kg/m3 to M, molec/cm3)
    ,KMAX_MID ,KMAX_BND, KCHEMTOP     ! Start and upper k for 1d fields
  use My_Aerosols_ml,       only : SEASALT
  use Landuse_ml,            only : water_fraction, ice_fraction
  use Par_ml,                only :  me& !!(me for tests)
                             ,gi0,gi1,gj0,gj1,IRUNBEG,JRUNBEG !hf VOL
  use PhysicalConstants_ml,  only :  AVOG, PI, GRAV
  use Radiation_ml,          only : PARfrac, Wm2_uE
  use Setup_1dfields_ml,     only : &
     xn_2d                &  ! concentration terms
    ,rcemis               &  ! emission terms
    ,rc_Rn222             &  ! for Pb210
    ,rcss                 &  !SeaS - sea salt
    ,rh, temp, tinv, itemp,pp      &  ! 
    ,amk, o2, n2, h2o           ! Air concentrations 
  use SeaSalt_ml,        only : SS_prod 
  use Tabulations_ml,    only :  tab_esat_Pa
  use TimeDate_ml,           only :  current_date, date
  use Volcanos_ml
  implicit none
  private
  !-----------------------------------------------------------------------!


  public :: setup_1d   ! Extracts results for i,j column from 3-D fields
  public :: setup_bio    ! Biogenic emissions
  public :: setup_rcemis ! Emissions  (formerly "poll")
  public :: reset_3d     ! Exports results for i,j column to 3-D fields

  real, dimension(NBVOC), public, save :: rcbio  !ispop and terpene

contains
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   subroutine setup_1d(i,j)

 !..   extracts data along one vertical column for input to chemical
 !     solver concentrations for chemistry......
 !
!DSGC ! Outputs, amk, o2k, rcairlig, ...

    integer, intent(in) :: i,j    ! coordinates of column

   !/* local

    integer           :: k, n, ispec, irc    ! loop variables
    real              :: qsat ! saturation water content


!DSGC    real ,dimension(KCHEMTOP:KMAX_MID) :: tinv, & ! Inverse of temp.
!DSGC        h2o, o2k     ! water, O2

    do k = KCHEMTOP, KMAX_MID
 
  !- MFAC - to scale from  density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

       amk(k) = roa(i,j,k,1) * MFAC  ! molecules air/cm3

       h2o(k) = max( 1.e-5*amk(k), &
                     q(i,j,k,1)*amk(k)*ATWAIR/18.0)

      ! nb. max function for h2o  used as semi-lagrangian scheme used
      ! in LAM50 (and HIRLAM) often gives negative H2O....   :-(

       pp(k) = A_mid(k) + B_mid(k)*ps(i,j,1)
       
       temp(k) = th(i,j,k,1)* Tpot_2_T( pp(k) )

       itemp(k) = nint( temp(k) -1.E-9)
! the "-1.E-9" is put in order to avoid possible different roundings on different machines. 

       qsat  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
       rh(k) = min( q(i,j,k,1)/qsat , 1.0) 

        ! 1)/ Short-lived species - no need to scale with M

         do n = 1, NSPEC_SHL
               xn_2d(n,k) = max(0.0,xn_shl(n,i,j,k))
         end do ! ispec

        ! 2)/ Advected species
        do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_2d(ispec,k) = max(0.0,xn_adv(n,i,j,k)*amk(k))
        end do ! ispec
 
        ! 3)/ Background species ( * CTM2 with units in mix. ratio)
        do n = 1, NSPEC_BGN
              xn_2d_bgn(n,k) = max(0.0,xn_bgn(n,i,j,k)*amk(k))
        end do ! ispec   

!DSGC ! setup weighting factor for hydrolysis  
!DSGC              f_Riemer(k)=96.*xn_2d(SO4,k)/( (96.*xn_2d(SO4,k))+(62.*xn_2d(aNO3,k)) )

   end do ! k

! Check that concentrations are not "contaminated" with NaN
   if ( .not.xn_2d(IXADV_NO2+NSPEC_SHL,KMAX_MID)+1>&
                        xn_2d(IXADV_NO2+NSPEC_SHL,KMAX_MID) ) then
      print *, "NANAN ", trim(species(IXADV_NO2+NSPEC_SHL)%name) 
   call CheckStop( "Detected non numerical concentrations (NaN)")
   end if

!DSGC o2k(:) = 0.21*amk(:)
   o2(:) = 0.21 *amk(:)
   n2(:) = amk(:) - o2(:)
!FAKE   o2(:) = 0.2095 *amk(:)
!FAKE   n2(:) = 0.7808 *amk(:)
   tinv(:) = 1./temp(:)


  ! 5 ) Rates  (!!!!!!!!!! NEEDS TO BE AFTER RH, XN, etc. !!!!!!!!!!)

!DSGC   rct(:,:)    = rcit(:,itemp(:))

   call set_rct_rates()

   !old: call set_rctroe_rates(tinv,amk,rctroe)

   call set_rcmisc_rates()
   !DSGC call set_rcmisc_rates(itemp,tinv,amk,o2k,h2o,rh,rcmisc)

  if ( DEBUG_SETUP_1DCHEM .and. debug_proc .and.  &
            i==debug_li .and. j==debug_lj .and. &
            current_date%seconds == 0 ) then
      write(*,"(a,f7.2,10es10.3)") " DEBUG_SETUP_1DCHEM ", &
            1.0/tinv(KMAX_MID), o2(KMAX_MID), &
            rcmisc(3,KMAX_MID), rcmisc(4,KMAX_MID), &
            rcmisc(10,KMAX_MID), rcmisc(11,KMAX_MID), &
            rcmisc(8,KMAX_MID), rcmisc(10,KMAX_MID)
      write(*,"(a,10es10.3)") " DEBUG_SETUP_1DCHEM RCT ", &
            rct(3,KMAX_MID), rct(4,KMAX_MID)
      write(*,"(a,10es10.3)") " DEBUG_SETUP_1DCHEM XN  ", &
        amk(KMAX_MID),  xn_2d(IXADV_O3+NSPEC_SHL,KMAX_MID), &
          xn_2d(IXADV_NO2+NSPEC_SHL,KMAX_MID)
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
     integer ::  iqrc,k,n, itot
     real    :: scaling, scaling_k
     real    :: eland   ! for Pb210  - emissions from land

    integer ::  i_help,j_help,i_l,j_l

! initilize
    rcemis(:,:)=0.    
    rcss(:,:) = 0.  !SeaS 
   
     do k=KEMISTOP,KMAX_MID

        do iqrc = 1, NRCEMIS
          itot = iqrc2itot( iqrc )
          rcemis(itot,k) = gridrcemis( iqrc ,k,i,j)
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
              !DSSRC rcemis(QRCVOL,k)=rcemis(QRCVOL,k)+rcemis_volc(volc_no)
              rcemis(SO2,k)=rcemis(SO2,k)+rcemis_volc(volc_no)
              !write(*,*)'Adding volc. emissions ' &
              !          ,rcemis_volc(volc_no),volc_no,&
              !          'to height=',k,'i,j',i_help,j_help
              !write(*,*)'TOT rcemis=',rcemis(QRCVOL,:)
           endif
        endif
     enddo
    
     endif ! VOLCANOES

    !/** lightning and aircraft ... Airial NOx emissions if required:

     if ( AIRNOX  ) then

        do k=KCHEMTOP, KEMISTOP,KMAX_MID

          rcemis(NO,k)  = rcemis(NO,k) &  
                              + 0.95 * (airn(k,i,j)+airlig(k,i,j))
          rcemis(NO2,k) = rcemis(NO2,k) &
                              + 0.05 * (airn(k,i,j)+airlig(k,i,j))

        enddo
        if ( DEBUG_SETUP_1DCHEM .and. debug_proc .and.  &
               i==debug_li .and. j==debug_lj ) write(*,"(a,10es10.3)") &
                 " DEBUG_SETUP_AIRNOX ", airn(KMAX_MID,i,j),airlig(KMAX_MID,i,j)

     end if ! AIRNOX

     !/** Add sea salt production
    
     if ( SEASALT  ) then

          do iqrc = 1, NSS
            rcss(iqrc,KMAX_MID) = SS_prod(iqrc,i,j)
          enddo

     endif

! Biogenics:
! Hard-coded for now, for isoprene only:
     rcemis(C5H8,KMAX_MID) =  rcemis(C5H8,KMAX_MID) + rcbio(BIO_ISOP)


     if ( FOREST_FIRES  .and. burning(i,j)  ) then

       call Fire_rcemis(i,j)

     endif  !ForestFires


!Mass Budget calculations
!   Adding up the emissions in each timestep


   scaling = dt_advec * xmd(i,j)* gridwidth_m*gridwidth_m / GRAV

   do k = KCHEMTOP,KMAX_MID

       
       scaling_k = scaling * (dA(k) + dB(k)*ps(i,j,1))/amk(k)

       do iqrc = 1, NSPEC_ADV
          
          itot = iqrc + NSPEC_SHL
          totem( iqrc ) = totem( iqrc ) + &
                rcemis( itot, k ) * scaling_k
          !if ( DEBUG_SETUP_1DCHEM .and. debug_proc .and.  &
          !    i==debug_li .and. j==debug_lj ) then
          ! write(6,"(a,2i3,es10.3,2i4)") "MASSEQV:", iqrc, rcemis( iqrc,k), qrc2ixadv(iqrc)
          !end if
       end do

   end do ! k loop 

  ! Soil Rn222 emissions from non-ice covered land, + water
  ! at rate of 1 atom/cm2/s

     eland = 1.0 - water_fraction(i,j) - ice_fraction(i,j)

! initialize, needed in My_Reactions
     rc_Rn222(:)=0.0     

! z_bnd is in m, not cm, so need to divide by 100.
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
  !  Output : rcbio added to rcemis - isoprene emissions for 1d column
  !
  !  Called from setup_ml, every  advection step.
  !----------------------------------------------------------------------------

!  input
  integer, intent(in) ::  i,j

!  local
  integer la,it2m,n,k,base,top,iclcat
  real clear

 ! Light effects added for isoprene emissions

  real            :: par   ! Photosynthetically active radiation
  real            :: cL    ! Factor for light effects
  real, parameter :: &
      CL1 = 1.066  , &    ! Guenther et al's params
      ALPHA = 0.0027      ! Guenther et al's params

  if ( NBVOC == 0  ) return   ! e.g. for ACID only



  it2m = nint(t2_nwp(i,j,1)-273.15-1.E-9)
! the "-1.E-9" is put in order to avoid possible different roundings on different machines.
  it2m = max(it2m,1)
  it2m = min(it2m,40)

  if(BIO_TERP > 0)  &
  rcbio(BIO_TERP) = emnat(i,j,BIO_TERP)*canopy_ecf(BIO_TERP,it2m)

  ! Isoprene has emissions in daytime only:
  rcbio(BIO_ISOP) = 0.0
  if ( Grid%izen <= 90) then

     ! Light effects from Guenther G93
      par = (Idirect(i,j) + Idiffuse(i,j)) * PARfrac * Wm2_uE
      cL = ALPHA * CL1 * par/ sqrt( 1 + ALPHA*ALPHA * par*par)

      rcbio(BIO_ISOP) = emnat(i,j,BIO_ISOP) &
             * canopy_ecf(BIO_ISOP,it2m) * cL
  endif
  if ( DEBUG_SETUP_1DBIO .and. debug_proc .and.  i==debug_li .and. j==debug_lj .and. &
         current_date%seconds == 0 ) then
     write(*,"(a5,2i4,4es12.3)") "DBIO ", current_date%day, &
      current_date%hour, par, cL, emnat(i,j,BIO_ISOP), rcbio(BIO_ISOP)
  end if

  end subroutine setup_bio

  !----------------------------------------------------------------------------
   subroutine reset_3d(i,j)

    integer, intent(in) :: i,j
    integer :: k, n, ispec, irc    ! loop variables
 
 
         do k = KCHEMTOP, KMAX_MID
 

!EGU - just testing units
!         xn2molem3 = 1.0/AVOG * 1.0e6
!         xn2ugC   = xn2molem3 * 12.0 * 1.0e6
!         xn2ugS   = xn2molem3 * 32.0 * 1.0e6
!         ugS2xn   = 1/xn2ugS = AVOG/(32.0*1.0e12)
!           xn_2d(SO4,k) =  0.5 * AVOG/(32.0*1.0e12)


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
