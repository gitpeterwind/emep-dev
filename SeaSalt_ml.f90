!_____________________________________________________________________________!
 module SeaSalt_ml
!SeaS

  !-----------------------------------------------------------------------! 
  ! Gets Inorganic Aerosols (amN, amS, SO4, ?NO3?) from MADE gas-phase
  ! chemistry, fraction additional masses to Aitken and accumulation modes
  ! based on the SO4 fraction. Particles number is unchanged.
  !st ? Split to SO4, NO3 and NH4? calculates new diameters  ??
  !-----------------------------------------------------------------------!
 use My_Emis_ml,           only : NSS, QSSFI, QSSCO
 use Par_ml,               only : MAXLIMAX,MAXLJMAX   ! => x, y dimensions
 use PhysicalConstants_ml, only : PI, AVOG
 use Par_ml,               only : me   ! for testing

  implicit none
  private

  public ::  SeaSalt_flux       

!  real, public, dimension(NSS,MAXLIMAX,MAXLJMAX) :: SS_prod
  real, public             :: SS_prod (NSS,MAXLIMAX,MAXLJMAX) = 0.
  real, public, save       :: n_to_mSS

  integer, parameter       :: ss_mod1 = 12, ss_mod2 = 3,   &  
                              nfin = 10, ncoa = 15,        &
                              SSdens = 2200.
  real, save                     :: param_SS(ss_mod2) = 0.  
  real, save, dimension(ss_mod2) :: radSS, dSS3
  real, save, dimension(ss_mod1) :: log_dp1, log_dp2 ,flux_help, dp3


  logical, private, save :: my_first_call = .true.

  contains

! <-------------------------------------------------------------------------------->

    subroutine SeaSalt_flux (i,j)

!!.. calculates generation of sea salt based on Maartinsson et al.(2003) 
!!.. for smaller particles and Monahan et al. for larger ones

 use Met_ml,               only : fm, roa, z_bnd, iclass
 use UKdep_ml,             only : landuse_ncodes, landuse_codes, landuse_data   !st sept,2004
 use DepVariables_ml,      only:  water                                         !st sept,2004
 use ModelConstants_ml,    only : KMAX_MID, KMAX_BND
 use GenSpec_tot_ml,       only : SSFI, SSCO
 use GenChemicals_ml,      only : species
 use PhysicalConstants_ml, only : CHARNOCK,GRAV,KARMAN, AVOG

 implicit none

   integer, intent(in) :: i,j    ! coordinates of column
   integer :: k, ii, jj, nlu, ilu, lu
   real    :: ustar, z0, z00, vind10, delz, zcoef, n2m, vind10_341,     &
              ss_flux(ss_mod1+ss_mod2), d3(ss_mod1+ss_mod2)  &
              ,prodM_ss(2),prodN_ss(2), water_frac
!//---------------------------------------------------

    if ( my_first_call ) then 

        call init_seasalt
    
        my_first_call = .false.

    end if !  my_first_call

    SS_prod(:,i,j) = 0.
    prodM_ss(:) = 0.
    prodN_ss(:) = 0.

   !ds ================
    if ( iclass(i,j) /= 0) return  !ds - faster to check here!
   !ds ================

!st sept,2004 - loop over the LU present in the grid
    nlu = landuse_ncodes(i,j)
      do ilu= 1, nlu
        lu      = landuse_codes(i,j,ilu)

!// only over water
!st sept,2004..  Obs!  All water is assumed here to be salt water
!                double checking with the old rough.170 data 

        !ds if ( water(lu) .and. iclass(i,j) == 0) then
        if ( water(lu) ) then
            water_frac = landuse_data (i,j,ilu)  ! grid fraction with water

          ustar = max(sqrt(fm(i,j,1)/roa(i,j,KMAX_MID,1)) , 1.e-5)

          z0 = CHARNOCK * ustar * ustar/GRAV
          z0 = max(z0,0.01)  ! u7.5vgc 
          vind10 = ustar/KARMAN * log(10./z0)

          delz = (z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND)) 
          zcoef = 1.e-6 / delz
          n2m = n_to_mSS * zcoef *AVOG / species(SSFI)%molwt *1.e-15
!pw         vind10341=vind10 ** (3.41)
          vind10_341=exp(log(vind10) * (3.41))

   do ii = 1, ss_mod1

!// sea salt flux [part/m2/s]

        ss_flux(ii) = flux_help(ii) * ( log_dp2(ii) - log_dp1(ii) )    &
!pw                                   * vind10 ** (3.41) 
                                   * vind10_341

        d3(ii) = dp3(ii)

   enddo

   do ii = 1, ss_mod2

!// sea salt flux [part/m2/s]

      jj = ii + ss_mod1

!pw        ss_flux(jj) = param_SS (ii)* vind10 **(3.41)
        ss_flux(jj) = param_SS (ii)*  vind10_341

!   if(bug) write(6,'(a15,i5,e13.4)') 'Flux Monah ->',jj, ss_flux(jj)

        d3(jj) = dSS3(ii)

   enddo
 
!// Sea salt particles Number production [ part/cm3/s]
!//.. ss_mod1= 12 (Martinsson et al.), ss_mod2= 3 (Monahan)

!..FINE ..

   do ii = 1 , nfin
      prodN_ss(QSSFI) = prodN_ss(QSSFI) + ss_flux(ii) *zcoef
      prodM_ss(QSSFI) = prodM_ss(QSSFI) + ss_flux(ii)* d3(ii) *zcoef

      SS_prod(QSSFI,i,j) = SS_prod(QSSFI,i,j)   &
                         + ss_flux(ii) * d3(ii) * n2m   &
                         * water_frac                     !st sept,2004 
   enddo

!..COARSE .. 
   do ii = nfin+1, ncoa
      prodN_ss(QSSCO) = prodN_ss(QSSCO) + ss_flux(ii) *zcoef
      prodM_ss(QSSCO) = prodM_ss(QSSCO) + ss_flux(ii)* d3(ii) *zcoef

      SS_prod(QSSCO,i,j) = SS_prod(QSSCO,i,j)   &
                         + ss_flux(ii) * d3(ii) * n2m   &
                         * water_frac                     !st sept,2004 
   enddo

! if(bug)  &
!  write(6,'(a20,2e15.4)') '>>Mass SS  >>', SS_prod(QSSFI,i,j), SS_prod(QSSCO,i,j)
! if(bug)  then
!  write(6,'(a20,2e15.4)') '>>N-SS prod >>', prodN_ss(QSSFI), prodN_ss(QSSCO)
!  write(6,'(a20,2e15.4)') '>>M-SS prod >>', prodM_ss(QSSFI), prodM_ss(QSSCO)
! endif                           
    endif  ! water_area > 0.
   enddo  ! LU classes

   end subroutine SeaSalt_flux

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine init_seasalt

  implicit none

  integer :: i, k
  real    :: a1, a2
  real, dimension(ss_mod2) :: range

!//=====
  real, parameter, dimension(5) ::    &
        c1 = (/-2.576e35,  5.932e28, -2.867e21, -3.003e13, -2.881e6 /),  &
        c2 = (/-2.452e33,  2.404e27, -8.148e20,  1.183e14, -6.743e6 /),  &
        c3 = (/ 1.085e29, -9.841e23,  3.132e18, -4.165e12,  2.181e6 /),  &

        d1 = (/ 7.188e37, -1.616e31,  6.791e23,  1.829e16,  7.609e8 /),  &
        d2 = (/ 7.368e35, -7.310e29,  2.528e23, -3.787e16,  2.279e9 /),  &
        d3 = (/-2.859e31,  2.601e26, -8.297e20,  1.105e15, -5.800e8 /)

!=== mikrometer
  real, parameter :: mkm  = 1.e-6,  mkm2 = 1.e-12 ,  &  
                     mkm3 = 1.e-18, mkm4 = 1.e-24  
!//.. for  Maartinsson, Nilsson .. parameterisation 
  real, parameter, dimension(ss_mod1) ::    &

  dp   = (/0.025,0.04, 0.075, 0.12,  0.17,  0.25, 0.36,  0.55,  0.7, 1.0,  1.6,  2.4  /), & ! dry diameter
  dp_1 = (/0.02, 0.03, 0.05,  0.10,  0.145, 0.20, 0.30,  0.419, 0.6, 0.8,  1.25, 2.0  /), & ! lower boundary
  dp_2 = (/0.03, 0.05, 0.10,  0.145, 0.20,  0.30, 0.419, 0.60,  0.8, 1.25, 2.0,  2.8  /)    ! upper boundary
!!dp_2 = (/0.06, 0.10, 0.20,  0.29 , 0.40,  0.60, 0.828, 1.20,  1.6, 2.5,  4.0,  5.6  /)    ! wet upper boundary


  real, dimension(ss_mod1) ::  dp2, dp4, a, b

!//... sea water temperature ... preliminary constant
  real, save :: Twater = 280.   

        log_dp1(:) = log10(dp_1(:))
        log_dp2(:) = log10(dp_2(:))

        dp2(:) = dp(:)  * dp(:)
        dp3(:) = dp2(:) * dp(:)
        dp4(:) = dp3(:) * dp(:)

!//..  for Monahan parameterisation

       ! radSS (1)= 6.3 /2.
       ! radSS (2)= 8.5 /2.
        radSS (1)= 6.0 /2.
        radSS (2)= 7.2 /2. 
        radSS (3)= 9.0 /2. 
!// cubic DRY diameter for SS particles (assumed SS_Dd=0.55*SS_Dw at Rh~90%)
        dSS3(1) =  ( 1.1*radSS(1)) **3
        dSS3(2) =  ( 1.1*radSS(2)) **3
        dSS3(3) =  ( 1.1*radSS(3)) **3

       ! range(1)= 1.4/2.    ! Monahan     5.6 - 7.0 wet diam
       ! range(2)= 3.0/2.    ! Monahan     7.0 - 10.0
        range(1)= 1.0/2.    ! Monahan     5.6 - 6.6 wet diam
        range(2)= 1.4/2.    ! Monahan     6.6 - 8.0
        range(3)= 2.0/2.    !             8.0 - 10.0

!//=====  Maartinsson, Nilsson .. parameterisation =======
   k = ss_mod1

  a(1:4) = c1(1)*dp4(1:4)*mkm4  + c1(2)*dp3(1:4)*mkm3  + c1(3)*dp2(1:4)*mkm2  &
         + c1(4)*dp(1:4)*mkm    + c1(5)
  a(5:7) = c2(1)*dp4(5:7)*mkm4  + c2(2)*dp3(5:7)*mkm3  + c2(3)*dp2(5:7)*mkm2  &
         + c2(4)*dp(5:7)*mkm    + c2(5)
  a(8:k) = c3(1)*dp4(8:k)*mkm4  + c3(2)*dp3(8:k)*mkm3  + c3(3)*dp2(8:k)*mkm2 &
         + c3(4)*dp(8:k)*mkm    + c3(5)

  b(1:4) = d1(1)*dp4(1:4)*mkm4  + d1(2)*dp3(1:4)*mkm3  + d1(3)*dp2(1:4)*mkm2  &
         + d1(4)*dp(1:4)*mkm    + d1(5)
  b(5:7) = d2(1)*dp4(5:7)*mkm4  + d2(2)*dp3(5:7)*mkm3  + d2(3)*dp2(5:7)*mkm2  &
         + d2(4)*dp(5:7)*mkm    + d2(5)
  b(8:k)= d3(1)*dp4(8:k)*mkm4   + d3(2)*dp3(8:k)*mkm3  + d3(3)*dp2(8:k)*mkm2 &
         + d3(4)*dp(8:k)*mkm    + d3(5)

  flux_help (:) = a(:) * Twater + b(:)

 if (me == 0) then
  write (6,*) '---   new SS -----'
  write (6,*)'     ====   A  ===='
  write (6,'(2(6e10.2/))')  (a (i), i=1,ss_mod1)
  write (6,*)'     ====   B  ===='
  write (6,'(2(6e10.2/))')  (b (i), i=1,ss_mod1)
  write (6,*)'     ====   F  ===='
  write (6,'(2(6e10.2/))')  (flux_help (i), i=1,ss_mod1)
 endif

  flux_help (:) =   flux_help (:) * 3.84e-6

!//======    Monahan   ===========================

   do i= 1, ss_mod2
         a1 = ( 0.380 - log10(radSS(i)) ) / 0.650
         a2 = 1.19 * exp(-a1*a1)

    param_SS(i) = 1.373 * radSS(i)**(-3) * range(i) *     &
                         (1.+ 0.057 * radSS(i)**1.05) * 10. **a2


 if (me == 0) then
  write(6,*)
  write(6,*) ' >> Sea Salt for <<', i,a1,a2,param_SS(i)
  write(6,*) '  TEST PARAMETERS *******', nfin, ncoa
 endif
   enddo

     n_to_mSS = PI * SSdens / 6.   ! SeaS
 
   end subroutine init_seasalt



 end module SeaSalt_ml

