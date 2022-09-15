  !-----------------------------------------------------------------------!
  !
  !     This is the main module that calculates J-values in the EMEP model,
  !     using the CloudJ modules stored in ModsCloudJ_mod.f90.
  !    
  !-----------------------------------------------------------------------!

MODULE CloudJ_mod
  
    use FJX_CMN_MOD  
    use FJX_SUB_MOD
    use FJX_INIT_MOD 
    use CLD_SUB_MOD,           only: CLOUD_JX

    use DefPhotolysis_mod      
    use GridValues_mod,        only: glat,glon,A_bnd,B_bnd,A_mid,B_mid,dA,dB
    use LandDefs_mod,          only: LandType, LandDefs
    use Landuse_mod,           only: LandCover
    use MetFields_mod,         only: ps,foundcloudwater,q,th,cc3d,  &
                                      roa, z_bnd, cw_met, &
                                      foundcloudicewater, ciw_met
    use Config_module,         only: KMAX_BND,KMAX_MID,KCHEMTOP,METSTEP,USES, &
                                      NPROC, IOU_INST, num_lev3d,lev3d, cloudjx_strat
    use Par_mod,               only: me,LIMAX,LJMAX
    use TimeDate_mod,          only: daynumber,current_date,date
    use ZchemData_mod,         only: rcphot, rcphotslice
    use Functions_mod,         only: Tpot_2_T
    use TimeDate_ExtraUtil_mod,only: date2string
    use SmallUtils_mod,        only: find_index
    use DerivedFields_mod,     only: d_3d, f_3d

    use Chemfields_mod,        only: xn_adv, xn_bgn, xn_shl, NSPEC_BGN, Dobson
    use ChemDims_mod,          only: NSPEC_ADV
    use ChemSpecs_mod
    use PhysicalConstants_mod, only: AVOG
                                  

    IMPLICIT NONE
    PUBLIC

    real, parameter ::  PI180 = 3.141592653589793d0/180.d0 
    real  GMTAU, ALBEDO, XGRD, YGRD, PSURF, SCALEH
    real  CF, PMID, ZDEL, ICWC, F1
    integer I,J,K,L,N
    integer NRAN, NJXX

!     --------------------params sent to CLOUD_JX-------------------------
    real                     :: U0,SZA,REFLB,SOLF,CLDCOR,FG0
    logical                    :: LPRTJ=.false.
    integer                    :: CLDFLAG,NRANDO,IRAN,LNRG,ICNT
    integer                    :: NICA,JCOUNT
    character*6, dimension(JVN_)  ::  TITLJXX

    integer, save :: photo_out_ix_no2_fj = -1
    integer, save :: photo_out_ix_o3a_fj = -1
    integer, save :: photo_out_ix_o3b_fj = -1
    
!---U0 = cos (SZA), SZA = solar zenith angle. Activate GIT?
!---REFLB = Lambertian reflectivity at the Lower Boundary
!---SOLF = solar flux factor for sun-earth distance
!---FG0 = scale for asymmetry factor to get equivalent isotropic (CLDFLAG=3 only)
!---LPRTJ = .true. = turn on internal print in both CLOUD_JX & PHOTO_JX
!--- P = edge press (hPa), Z = edge alt (m), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!--- R = layer rel.hum.(fraction)
!---LWP/IWP = Liquid/Ice water path (g/m2)
!---REFFL/REFFI = R-effective(microns) in liquid/ice cloud
!---CLF = cloud fraction (0.0 to 1.0)
!---CLDIW = integer denoting cloud in layer: 1=water, 2=ice, 3=both
!---AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
!---  aerosols are dimensioned with up to AN_ different types in an ICA layer
!---L1_ = parameter, dim of profile variables, L_+1 for top (non CTM) layer
!---AN_ = parameter, dim of number of aerosols being passed (25 standard)
!---VALJXX = J-values from CLOUD_JX & PHOTO_JX
!---JVN_ = dim of max number of J-s reported out (in the order of fast-JX, not CTM)
!---CLDFLAG = integer index for type of cloud overlap
!---CLOUD_JX:   different cloud schemes
!---CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
!       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
!---NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
!---IRAN = starting index for random number selection
!---LNRG = flag for setting max-ran overlap groups:
!---     =0   break max overlap groups at cloud fraction = 0
!---     =3   else fixed 3 layers (1:9, 9:last LWcloud, LWclud+1:LTOP)
!---     else(=6) fixed correlated length max-overlap layers
!---NICA = total number of ICAs

!-----these are the key results coming out of fast-JX core
!-----they need to be dimensioned here and passed to fast-JX to fill.

CONTAINS

SUBROUTINE setup_phot_cloudj(i_emep,j_emep,errcode,mode)

    !calculate rcphot for one ICA in the emep model
    IMPLICIT NONE

    integer, intent(in) :: i_emep, j_emep
    integer, intent(inout) :: errcode
    integer,  intent(in) :: mode
    integer :: ilu, lu, nr_local
    real :: Pres_mid, temperature, swp

    integer IDAY, JLAT, ILON

    real, allocatable, save :: CWIC(:), CWLC(:)
    real, allocatable, save :: DeltZ(:),CLWP(:),CIWP(:)
    real, allocatable, save :: AER1(:),AER2(:),AER3(:),AER4(:),AER5(:)
    integer,allocatable, save :: NAA1(:), NAA2(:), NAA3(:), NAA4(:), NAA5(:)
    real, allocatable, save :: DUST_F(:),DUST_C(:),SULF(:),SEAS_F(:),SEAS_C(:)

    logical, save :: first_call=.true.

    ! ICA input to CloudJ. NB: L_in sets the number of levels read in by CLOUD_JX
    real, allocatable, save :: PPP(:), ZZZ(:), TTT(:), DDD(:), RRR(:)
    real, allocatable, save :: REFFI(:), CLF(:), OOO(:), LWP(:), REFFL(:), IWP(:)
    real, allocatable, save :: AERSP(:,:)

    integer, allocatable, save :: CLDIW(:)
    integer, allocatable, save :: NDXAER(:,:)
    real,  allocatable, save :: VALJXX(:,:)

    real :: col_o3

    ! declarations used in reading/handling satellite stratospheric ozone data
    integer :: la, lo, w, z
    integer, save :: dim_lon, dim_lat, dim_alt, OZ_TOP, oz_month=-999, IO_OZONE=78
    integer, save :: year, month, naerosol=5 ! naerosol sets number of included EMEP aerosol, must be set 
    real, allocatable, save :: lon_ozone(:), lat_ozone(:), alt_ozone(:)
    real, allocatable, save :: temp_ozone(:,:,:), pres_ozone(:,:,:), ozon_ozone(:,:,:)
    real, allocatable, save :: temp_obs(:,:,:), pres_obs(:,:,:), ozon_obs(:,:,:)
    character(len=150), save  :: fname_ozone
    real,  save :: dlat, dlon, dalt
    
    !---fast-JX:  INIT_JX is called only once to read in & store all fast-JX data: 
    !             also sets up random sequence for cloud-JX. 
    ! JVN_ = 101 (max number of J values read in from CloudJ)
    ! NJXX = number of derived Jvalues, set in RD_XXX 

    ! if(me==0.and.i_emep==1.and.j_emep==1)write(*,*) 'this is a test', LWEPAR 

    ! set time-variables; nr_local = 1 uses current meteorology timestep 
    nr_local=1
    if(mode>0)nr_local=mode
    iday = daynumber  
    year = current_date%year
    month = current_date%month
    gmtau = current_date%hour + current_date%seconds/3600.0
    if(nr_local==2)gmtau = gmtau + METSTEP ! use next met tstep

    ! read monthly lat-lon-altitude ozone satellite observations
    if(oz_month/=month) then 
          ! ##########################################################################################
          ! Ozone satellite observation files are preprocessed using the script Sat_Ozone_datsaver.py.
          ! Data is defined for altitude levels strictly above 100 hPa, so as not to                  
          ! conflict with the model top of the (current) EMEP model.                                  
          !                                                                                                            
          ! Data is based on measurements from six satellite missions between yrs 2005-2021:   
          !     https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-ozone-v1?tab=form
          !
          ! ##########################################################################################

          ! define path based on model year & month and read in monthly lon-lat-alt satellite obs. data        
          if(year < 2005 .or. year > 2021) then 
                fname_ozone = trim(cloudjx_strat)//trim('clim_')//date2string("MM",current_date)//trim('.dat')
          else
                fname_ozone = trim(cloudjx_strat)//date2string("YYYYMM",current_date)//trim('.dat')
          endif
          if(me==0)write(*,*) 'Opening satellite O3/T obs. file: ', fname_ozone
          
          open(IO_OZONE, file=fname_ozone,form='unformatted')
          read(IO_OZONE) dim_lon, dim_lat, dim_alt

          if(first_call) then 
                ! i/o arrays necessary for strat ozone files
                allocate(lon_ozone(dim_lon))
                allocate(lat_ozone(dim_lat))
                allocate(alt_ozone(dim_alt))

                allocate(temp_ozone(dim_lon,dim_lat,dim_alt))
                allocate(pres_ozone(dim_lon,dim_lat,dim_alt))
                allocate(ozon_ozone(dim_lon,dim_lat,dim_alt))

                allocate(temp_obs(LIMAX,LJMAX,dim_alt))
                allocate(pres_obs(LIMAX,LJMAX,dim_alt))
                allocate(ozon_obs(LIMAX,LJMAX,dim_alt))

                ! total number of levels with O3/T data
                OZ_TOP = KMAX_MID + dim_alt 

                ! i/o arrays necessary for cloudj
                allocate(PPP(OZ_TOP+1),ZZZ(OZ_TOP+1)) 
                allocate(TTT(OZ_TOP),DDD(OZ_TOP),RRR(OZ_TOP),CLF(OZ_TOP))
                allocate(OOO(OZ_TOP),LWP(OZ_TOP),IWP(OZ_TOP))
                allocate(REFFI(OZ_TOP),CLDIW(OZ_TOP),REFFL(OZ_TOP))
                allocate(NDXAER(OZ_TOP,naerosol),AERSP(OZ_TOP,naerosol))
                allocate(VALJXX(OZ_TOP-1, JVN_))

                allocate(CWIC(OZ_TOP),CWLC(OZ_TOP),DeltZ(OZ_TOP),CLWP(OZ_TOP))
                allocate(CIWP(OZ_TOP),AER1(OZ_TOP),AER2(OZ_TOP),AER3(OZ_TOP))
                allocate(AER4(OZ_TOP),AER5(OZ_TOP),DUST_F(OZ_TOP),DUST_C(OZ_TOP))
                allocate(SULF(OZ_TOP),SEAS_F(OZ_TOP),SEAS_C(OZ_TOP),NAA1(OZ_TOP))
                allocate(NAA2(OZ_TOP),NAA3(OZ_TOP),NAA4(OZ_TOP),NAA5(OZ_TOP))
          
                if(USES%EtaCOORDINATES.and.A_bnd(1) < 1.e4) write(*,*) 'Warning: CloudJ stratospheric O3/T UNDEFINED! ', &
                      'Top model level above 100 hPa.'
          end if

          read(IO_OZONE) lon_ozone ! -180 to 180 degrees, as in EMEP
          read(IO_OZONE) lat_ozone ! -90 to 90 degrees, as in EMEP
          read(IO_OZONE) alt_ozone ! altitude levels in km

          dlon = lon_ozone(2) - lon_ozone(1) ! grid spacing in degrees
          dlat = lat_ozone(2) - lat_ozone(1) 
          dalt = alt_ozone(2) - alt_ozone(1) ! vertical spacing in km

          read(IO_OZONE) temp_ozone ! temperature (K)
          read(IO_OZONE) pres_ozone ! pressure (hPa)
          read(IO_OZONE) ozon_ozone ! ozone (molec/cm2)

          close(IO_OZONE)     

          do w=1,LIMAX ! map lat(LJMAX) and lon(LIMAX) indices and populate the obs arrays
                do z=1,LJMAX 
                      lo = max(1,int( 1 / dlon * glon(w,z) + 179.999 / dlon + 1))
                      la = max(1,int( 1 / dlat * glat(w,z) + 89.999 / dlat + 1))

                      lo = min(dim_lon, lo) ! avoid any overshoot (e.g. when lon = 180.25 deg)
                      la = min(dim_lat, la)

                      temp_obs(w,z,:) = temp_ozone(lo,la,:)
                      pres_obs(w,z,:) = pres_ozone(lo,la,:)
                      ozon_obs(w,z,:) = ozon_ozone(lo,la,:)
                end do
          end do
    end if ! oz_month /= month
    
    oz_month = month

    ! initialize requested (in config file) J-value output arrays
    if(first_call)then 
          call INIT_FJX(TITLJXX,JVN_,NJXX) 

          if(allocated(f_3d))then ! if output fields requested in config file, set index
            photo_out_ix_no2_fj = find_index("D3_J(NO2)", f_3d(:)%subclass)
            photo_out_ix_o3a_fj = find_index("D3_J(O3a)", f_3d(:)%subclass)
            photo_out_ix_o3b_fj = find_index("D3_J(O3b)", f_3d(:)%subclass)
          endif

          if(photo_out_ix_no2_fj>0.and.me==0)write(*,*)'will output jNO2', photo_out_ix_no2_fj
          if(photo_out_ix_o3a_fj>0.and.me==0)write(*,*)'will output jO3a', photo_out_ix_o3a_fj
          if(photo_out_ix_o3b_fj>0.and.me==0)write(*,*)'will output jO3b', photo_out_ix_o3b_fj

          write(*,*)'Found cloud water (1) and ice (2): ', foundcloudwater, foundcloudicewater
    endif

    !inputs for SOLAR_JX 
    YGRD = glat(i_emep,j_emep)*PI180
    XGRD = glon(i_emep,j_emep)*PI180
    PSURF = ps(i_emep,j_emep,nr_local)/100.0 !Pa-> hPa

    ILON = int(XGRD/PI180)
    JLAT = int(YGRD/PI180)
    IRAN = 13+ILON+3*JLAT+5*(YEAR-1900)+7*IDAY + 11*nint(GMTAU) ! random number based on year/day/hour

    ! initiliaze cloud fraction and water path columns
    CLF(:)  = 0.
    CLWP(:) = 0.
    CIWP(:) = 0.
    CWIC(:) = 0.
    CWLC(:) = 0.
    RRR(:) = 0.

    AER1(:) = 0.
    AER2(:) = 0.
    AER3(:) = 0.
    AER4(:) = 0.
    AER5(:) = 0.

    NAA1(:) = 0
    NAA2(:) = 0
    NAA3(:) = 0
    NAA4(:) = 0
    NAA5(:) = 0

    !use fastj vertical direction, i.e. L largest at top, 1 at surface 
    L=0

    ! EMEP mid-levels 
    do k=kmax_mid,1,-1
          L = L+1
          Pres_mid = A_mid(k) + B_mid(k)*ps(i_emep,j_emep,nr_local)
    
          !potential -> absolute temperature
          temperature = th(i_emep,j_emep,k,nr_local)* Tpot_2_T( Pres_mid )
          TTT(L) = temperature

          !specific -> relative humidity
          swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65)) !saturation water pressure
          RRR(L) = q(i_emep,j_emep,k,nr_local)*(Pres_mid)/0.622/swp

          !cloud cover stored as fraction from 0 to 1 (line 652 Met_mod.f90)
          DeltZ(L) = z_bnd(i_emep,j_emep,k)-z_bnd(i_emep,j_emep,k+1)
          CF = cc3d(i_emep,j_emep,k,1)

          !if clouds are present, populate cloud liquid and ice water content arrays
          if(CF > 0.00001d0) then
                CLF(L) = CF
                if(foundcloudwater) then
                   CWLC(L) = cw_met(i_emep,j_emep,k,1) / CLF(L)  ! kg/kg per cloud area 
                   CLWP(L) = CWLC(L) * 1000. * roa(i_emep,j_emep,k,1) * DeltZ(L)  ! water path [g/m2]
                endif  
                if(foundcloudicewater) then
                   CWIC(L) = ciw_met(i_emep,j_emep,k,1) / CLF(L) ! kg/kg per cloud area 
                   CIWP(L) = CWIC(L) * 1000. * roa(i_emep,j_emep,k,1) * DeltZ(L)  ! water path [g/m2]
                endif
          end if

          ! for aerosols, positive indices set std atmos defs. 
          
          ! name(a12)_| R-eff   rho| notes
          ! 1 Str-Bkgd| 0.221 1.630| backgrd strat sulf
          ! 2 Str-Volc| 0.386 1.630| volcanic strat sulf
          ! 3 UT-sulf1| 0.166 1.769| upper-trop sulf1
          ! 4 UT-sulf2| 0.166 1.769| upper-trop sulf2
          ! 5 UT-sulfM| 0.140 1.769| upper-trop sulf (UMich)
          ! 6 UM-BC1  | 0.140 1.500| UMich w/Mie code
          ! 7 UM-BC2  | 0.140 1.500| UMich w/Mie code
          ! 8 UM-BB08C| 0.149 1.230| UMich w/Mie code 8%BC 0%RH
          ! 9 UM-FF04C| 0.140 1.212| UMich w/Mie code 4%BC 0%RH
          !10 UM-FF10C| 0.140 1.230| UMich w/Mie code 10%BC 0%RH
          !11 HDust.15| 0.150 2.600| Harvard R.V.Martin generated (original)
          !12 HDust.25| 0.250 2.600
          !13 HDust.40| 0.400 2.600
          !14 HDust.80| 0.800 2.600
          !15 HDust1.5| 1.500 2.600
          !16 HDust2.5| 2.500 2.600
          !17 HDust4.0| 4.000 2.600

          ! negative indices give UMich aerosol definitions, which include RH effects

          ! 1 SULF    2 SS-1    3 SS-2    4 SS-3    5 SS-4    6 DD-1    7 DD-2
          ! 8 DD-3    9 DD-4    10 FF00   11 FF02   12 FF04   13 FF06   14 FF08
          ! 15 FF10   16 FF12   17 FF14   18 BB00   19 BB02   20 BB04   21 BB06
          ! 22 BB08   23 BB10   24 BB12   25 BB14   26 BB16   27 BB18   28 BB20
          ! 29 BB22   30 BB24   31 BB26   32 BB28   33 BB30

          !species            mol wts
          !OD                 16.0000
          !IXADV_Dust_road_f 200.0000
          !IXADV_Dust_wb_f   200.0000
          !IXADV_Dust_sah_f  200.0000
          !IXADV_Dust_road_c 200.0000
          !IXADV_Dust_wb_c   200.0000
          !IXADV_Dust_sah_c  200.0000

          !IXADV_pSO4f          96.0000  
          !IXADV_pSO4c          96.0000
          !IXADV_SeaSalt_f      58.0000  
          !IXADV_SeaSalt_c      58.0000  
          !IXADV_ffire_BC       12.0000  
          !IXADV_ffire_c        12.0000
          !IXADV_Ash_f          12.0000 
          !IXADV_Ash_c          12.0000
          !IXADV_POM_f_wood     20.4000
          !IXADV_POM_c_wood     20.4000
          !IXADV_POM_f_ffuel    15.0000
          !IXADV_POM_c_ffuel    15.0000
          !IXADV_NO3_f          62.0000  
          !IXADV_NO3_c          62.0000 
          !IXADV_NH4_f          18.0000

          !density from kg to gram, divided by mean molar mass, times altitude, times molar mass dust
          DUST_F(L) = (xn_adv(IXADV_Dust_road_f,i_emep,j_emep,k) &
                    +  xn_adv(IXADV_Dust_wb_f,i_emep,j_emep,k)   &
                    +  xn_adv(IXADV_Dust_sah_f,i_emep,j_emep,k)) &
                    * 1000 * roa(i_emep,j_emep,k,1) / 28.971600  &
                    * DeltZ(L) * 200 

          DUST_C(L) = (xn_adv(IXADV_Dust_road_c,i_emep,j_emep,k) &
                    +  xn_adv(IXADV_Dust_wb_c,i_emep,j_emep,k)   &
                    +  xn_adv(IXADV_Dust_sah_c,i_emep,j_emep,k)) &
                    * 1000 * roa(i_emep,j_emep,k,1) / 28.971600  &
                    * DeltZ(L) * 200

          sulf(l) = 0.

          ! SULF(L)   = (xn_adv(IXADV_pSO4f,i_emep,j_emep,k)  & not all chemistry has sulf
          !           +  xn_adv(IXADV_pSO4c,i_emep,j_emep,k)) &
          !           * 1000 * roa(i_emep,j_emep,k,1) / 28.971600  &
          !           * DeltZ(L) * 96

          SEAS_F(L) = xn_adv(SeaSalt_f,i_emep,j_emep,k) &
                     * 1000 * roa(i_emep,j_emep,k,1) / 28.971600  &
                     * DeltZ(L) * 58

          SEAS_C(L) = xn_adv(SeaSalt_c,i_emep,j_emep,k) &
                     * 1000 * roa(i_emep,j_emep,k,1) / 28.971600  &
                     * DeltZ(L) * 58
          
          ! naerosol variable MUST be set with the number of used aeorsol species 
          AER1(L)=DUST_F(L) ! aerosol [g/m2]
          NAA1(L)=-8        ! aerosol index

          AER2(L)=DUST_C(L) ! aerosol [g/m2]
          NAA2(L)=-9        ! aerosol index

          AER3(L)=SULF(L)   ! aerosol [g/m2]
          NAA3(L)=-1        ! aerosol index

          AER4(L)=SEAS_F(L) ! aerosol [g/m2]
          NAA4(L)=-4        ! aerosol index

          AER5(L)=SEAS_C(L) ! aerosol [g/m2]
          NAA5(L)=-5        ! aerosol index   
    end do ! EMEP temperature, relative humidity, cloud water path & aerosol

    L=0
    do k=KMAX_BND,1,-1 ! EMEP level interface pressure
          L=L+1
          PPP(L) =A_bnd(k)/100.0 + B_bnd(k)*PSURF ! PSURF in hPa
    end do
    
    L=0 ! EMEP ozone, derived using interface pressure
    do k=KMAX_MID,1,-1
          L=L+1
          OOO(L) = xn_adv(IXADV_O3,i_emep,j_emep,k)*(PPP(L)-PPP(L+1))*MASFAC
    end do

!------------------------------------------------------------------------------------
! Set layers above the EMEP model
!------------------------------------------------------------------------------------
   
    K = 0 ! populate stratospheric ozone and temperature
    do L=KMAX_BND,OZ_TOP
          K = K + 1
          OOO(L) = ozon_obs(i_emep,j_emep,K)
          TTT(L) = temp_obs(i_emep,j_emep,K)
          PPP(L+1) = pres_obs(i_emep,j_emep,K) 
    end do
    
    !-----------------------------------------------------------------------
    !---fast-JX:  SOLAR_JX is called only once per grid-square to set U0, etc. 
    call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD,SZA,U0,SOLF)

    LWP(:)  = 0.d0       ! liquid water path (g/m2)
    IWP(:)  = 0.d0       ! ice water path (g/m2)
    AERSP(:,:)  = 0.d0   ! aerosol path (g/m2)
    NDXAER(:,:) = 0      ! aerosol index type
    CLDIW(:) = 0

    REFFL(:) = 0.d0      ! R-effective(microns) in liquid cloud
    REFFI(:) = 0.d0      ! R-effective(microns) in ice cloud

    ZZZ(1) = 16.d5*log10(1013.25d0/PPP(1))
    
    Dobson(i_emep,j_emep) = 0. ! this array may now also be used simply to test
    col_o3=0.

    do L = 1,OZ_TOP 
          DDD(L) = (PPP(L)-PPP(L+1))*MASFAC
          SCALEH = 1.3806d-19*MASFAC*TTT(L)
          ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH ) ! can be changed to use DeltZ, to do
          
          ! CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
          if (CWLC(L) > 1.d-11) CLDIW(L) = 1
          if (CWIC(L) > 1.d-11) CLDIW(L) = CLDIW(L) + 2

          if (CWLC(L) > 1.d-12) then            ! [kg/kg]
                LWP(L) = CLWP(L)                   ! [g/m2]
                PMID = 0.5d0*(PPP(L)+PPP(L+1))
                F1   = 0.005d0 * (PMID - 610.d0)
                F1   = min(1.d0, max(0.d0, F1))
                REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)    
          endif

          if (CWIC(L) > 1.d-12) then            ! [kg/kg]
                IWP(L) = CIWP(L)                   ! [g/m2]
                ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! [m]
                ICWC = IWP(L)/ZDEL                 ! [g/m3]
                REFFI(L) = 164.d0 * (ICWC**0.23d0)
          endif

          NDXAER(L,1) = NAA1(L) ! populate aerosol arrays
          AERSP(L,1)  = AER1(L)
          NDXAER(L,2) = NAA2(L)
          AERSP(L,2)  = AER2(L)
          NDXAER(L,3) = NAA3(L)
          AERSP(L,3)  = AER3(L)
          NDXAER(L,4) = NAA4(L)
          AERSP(L,4)  = AER4(L)
          NDXAER(L,5) = NAA5(L)
          AERSP(L,5)  = AER5(L)
          
          if (L > KMAX_MID) &  ! ozone column above EMEP model layers
                col_o3 = col_o3 + OOO(L)

          dobson(i_emep,j_emep) = dobson(i_emep,j_emep) + OOO(L)/2.687e16 !1 DU = 2.687e16 molecules of O3 per square centimetre
    end do

    ! albedo read in as %. Get land-name, then read albedo and convert to fraction 
    ALBEDO = 0.0
    do ilu=1, LandCover(i_emep,j_emep)%ncodes
          lu = LandCover(i_emep,j_emep)%codes(ilu)
          ALBEDO = ALBEDO + LandDefs(lu)%Albedo*0.01*LandCover(i_emep,j_emep)%fraction(ilu)
    end do
    REFLB = ALBEDO
                      
    !--CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
    !       CLDFLAG = 1  :  Clear sky J's
    !       CLDFLAG = 2  :  Averaged cloud cover
    !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
    !       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
    !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
    !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
    !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin) ! recommended Prather 2015
    !       CLDFLAG = 8  :  Calcluate J's for ALL ICAs (up to 20,000 per cell!)
    CLDFLAG = 3        
    
    ! asymmetry factor equivalent isotropy: CLDFLAG = 3 only   
    FG0     = 1.1

    ! for CLDFLAG > 3 only
    CLDCOR  = 0.33
    NRANDO  = 5
    LNRG    = 6

    LPRTJ = .false.
    if(me==0.and.first_call) LPRTJ = .false. ! may want to choose verbose on first call 
    
    !=======================================================================
    ! outputs Jvalues (VALJXX), NICA and JCOUNT
    !
    ! inputs PPP in Pascal, ZZZ in cm, TTT in Kelvin, DDD in molec/cm2,
    ! RRR as RH (0-100%), OOO in molec/cm2, LWP in g/m2, IWP in g/m2
    !=======================================================================
    call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT,             &
                DDD,RRR,OOO,LWP,IWP,REFFL,REFFI,CLF,CLDCOR,CLDIW, &
                AERSP,NDXAER,OZ_TOP,naerosol,VALJXX,JVN_,                    &
                CLDFLAG,NRANDO,IRAN,LNRG,NICA,JCOUNT)
    !=======================================================================
    
    !---map multiplication factors to J-values from fast-JX using JIND & JFACTA.
    !---For example, JFACTA(1) from photochemistry scheme maps to FJX(4)

    ! apply multiplication factors for reactions that are based on other reactions
    do L = 1,OZ_TOP-1
          do J = 1,NRATJ
                ! only if J-value defined (memory issues otherwise)
                if(JIND(J)>0) VALJXX(L,JIND(J)) = VALJXX(L,JIND(J))*JFACTA(J)
          end do
    end do
    
    ! populate 3D Jvalue array with CloudJ values for the photochemical reactions in EMEP
    if(.not.(allocated(rcphotslice))) allocate(rcphotslice(NRCPHOTextended,KCHEMTOP:KMAX_MID,LIMAX,LJMAX))
    
    do L=1,KMAX_BND-KCHEMTOP ! KMAX_BND = 21, KMAX_MID = 20, KCHEMTOP = 2
          ! reactions having a 1-to-1 correspondence with tabulated reactions. Note: All of these CloudJ reactions have been checked to be non-zero.
          rcphotslice(IDAO3,kmax_bnd-L,i_emep,j_emep)    = VALJXX(L,JIND(3))   ! JIND(J))  *JFACTA(J)   ZPJ(L,3)!3 O3   PHOTON    O2    O(total)   1.000 /O3    /
          rcphotslice(IDBO3,kmax_bnd-L,i_emep,j_emep)    = VALJXX(L,JIND(4))   ! ZPJ(L,4)!4  O3        PHOTON    O2        O(1D)                   1.000 /O3(1D) /           
          rcphotslice(IDNO2,kmax_bnd-L,i_emep,j_emep)    = VALJXX(L,JIND(9))   !  9 NO2       PHOTON    N2        O                       1.000 /NO2   /
          rcphotslice(IDH2O2,kmax_bnd-L,i_emep,j_emep)   = VALJXX(L,JIND(7))   !  7 H2O2      PHOTON    OH        OH                      1.000 /H2O2  /
          rcphotslice(IDHNO3,kmax_bnd-L,i_emep,j_emep)   = VALJXX(L,JIND(15))  ! 15 HNO3      PHOTON    NO2       OH                      1.000 /HNO3  /
          rcphotslice(IDACH2O,kmax_bnd-L,i_emep,j_emep)  = VALJXX(L,JIND(5))   !  5 H2CO      PHOTON    HCO       H                       1.000 /H2COa /
          rcphotslice(IDBCH2O,kmax_bnd-L,i_emep,j_emep)  = VALJXX(L,JIND(6))   !  6 H2CO      PHOTON    CO        H2                      1.000 /H2COb /
          rcphotslice(IDHONO,kmax_bnd-L,i_emep,j_emep)   = VALJXX(L,JIND(14)) ! 14 HNO2      PHOTON    NO       OH    
          rcphotslice(IDHO2NO2,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(16)) / 0.667 ! emep ratio is 0.667 for NO2 + HO2 
          rcphotslice(IDNO3,kmax_bnd-L,i_emep,j_emep)    = VALJXX(L,JIND(11)) + VALJXX(L,JIND(12)) ! emep ratio is 0.873 to 0.127, fastj ratio is 0.886 to 0.114
          rcphotslice(IDCH3O2H,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(8))  ! Methyl hydroperoxide. Used in many photolysis reactions. See Blitz et al. 2005

          ! pending:
          ! CH3COCHO paper for reference, but IDRCOHCO currently not in use, only through IDRCOCHO = IDHCOHCO = 11. Can be enabled separately by defining in DefPhot.f9
          rcphotslice(IDRCOHCO,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(64))!64 CH3COCHO  PHOTON    CH3CO     CO                      1.000 /MGlyxl/

          !emep CH3CHO -> CH3O2+HO2+CO -> CH3 + O2 + HO2 + CO, fastj: CH3CHO -> CH3 + HCO -> CH3 + HCO2 (+ O2) + CO; Moortgat et al. 2010
          rcphotslice(IDCH3CHO,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(54)) !54 CH3CHO    PHOTON    CH3       HCO                     1.000 /ActAld/


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! not the same, disable?
          ! IDCH3COX = IDMEK. Tabulated is about a factor of 7 larger. 
          rcphotslice(IDCH3COX,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(61))+VALJXX(L,JIND(62))

          ! IDCH3COY = C4H6O2 groups? Bi-acetone?
          rcphotslice(IDCH3COY,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(69))!69 CH3COCH3  PHOTON    CH3       CH3       CO            1.000 /Acet-b/ ?

          !emep: IDCHOCHO = IDHCOHCO, GLYOX=HCOHCO -> CO(1.9) + HO2(0.5) + HCHO(0.1). Not the same reactions, it seems! 
          rcphotslice(IDHCOHCO,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(66))!66 CHOCHO    PHOTON    H2        CO        CO            1.000 /Glyxlb/


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! not in use in current photochemistry scheme
          rcphotslice(IDACETON,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(68)) !not used !68 CH3COCH3  PHOTON    CH3CO     CH3                     1.000 /Acet-a/ 
          !note N2O5: generally low concentrations during daytime anyway due to NO3 photolysis (Wayne et al 2011)         
          rcphotslice(IDN2O5,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(13))   !not used !13 N2O5      PHOTON    NO2       NO3                     1.000 /N2O5  /


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! to add: PAN ?
    end do

    first_call=.false.
    LPRTJ=.false.

end subroutine setup_phot_cloudj


subroutine write_jvals(i_emep,j_emep)

    ! this routine writes the desired 3D J-value output. Indices are set either by 
    ! setp_phot_cloudj (CloudJ) or setup_phot (DefPhotolysisMod.f90, tabulated)

    integer, intent(in) :: i_emep,j_emep
    
    if(photo_out_ix_no2>0)then
          d_3d(photo_out_ix_no2,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDNO2,max(KCHEMTOP,lev3d(1:num_lev3d))) !WARNING: rcphot defined only up to KCHEMTOP!
    endif

    if(photo_out_ix_o3a>0)then
          d_3d(photo_out_ix_o3a,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDAO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_o3b>0)then
          d_3d(photo_out_ix_o3b,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDBO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_no2_fj>0)then
          d_3d(photo_out_ix_no2_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDNO2,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_o3a_fj>0)then
          d_3d(photo_out_ix_o3a_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDAO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_o3b_fj>0)then
          d_3d(photo_out_ix_o3b_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDBO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

end subroutine write_jvals


endmodule CloudJ_mod














