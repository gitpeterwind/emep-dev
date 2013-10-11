! <AOD_PM_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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


module AOD_PM_ml
!-----------------------------------------------------------------------! 
! Calculates Aerosol Optical Depth (AOD) for 0.5 um radiation based on 
! aerosol mass concentrations and specific extinction cross-sections 
! based on Tegen et al. JGR (1997) and Kinne et al., ACP (2005)
! (implicit assumption on wetted aerosols)
!-----------------------------------------------------------------------!
use Chemfields_ml,        only: AOD, Extin_coeff
use ChemChemicals_ml,     only: species  
use ChemGroups_ml,        only: AOD_GROUP
use GridValues_ml,        only: i_fdom, j_fdom
use MetFields_ml,         only: z_bnd
use ModelConstants_ml,    only: KMAX_MID, KCHEMTOP
use Par_ml,               only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use PhysicalConstants_ml, only: AVOG
use Setup_1dfields_ml,    only: xn_2d, rh
use SmallUtils_ml,        only: find_index

implicit none
private
!-----------------------------------------------------------------------!
!// Subroutines
public :: AOD_calc, AOD_Ext
!// Functions
public :: wetExtC

real, parameter :: lambda = 0.55e-6
contains
! <---------------------------------------------------------->
function wetExtC(ntot,gtot,rh,debug) result(extC)
!.. Calculate Rh dependent specific extinction (or mass extinction efficiency)
!..   according to Chin et.al (J. Atm.Sci., 59, 2001)
  integer, intent(in) :: ntot
  integer, dimension(ntot), intent(in) :: gtot
  real, intent(in) :: rh
  logical, intent(in) :: debug  
  real, dimension(ntot) :: extC
  real, parameter :: &
    rhoSO4=1.6, rhoOC=1.8, rhoEC=1.0, rhoDU=2.6, rhoSS=2.2 , &
    Reff_SO4 = 0.156, Reff_OC = 0.087, Reff_EC = 0.039, &
    Reff_DUf = 0.80,  Reff_DUc = 4.5, Reff_SSf = 0.80, Reff_SSc = 5.73
  integer, parameter :: NumRH=7
  real, parameter, dimension(NumRH) ::    &
    RelHum = (/ 0.0,  0.5, 0.7, 0.8, 0.9, 0.95, 0.99 /),  &
    GF_SO4 = (/ 1.0,  1.4, 1.5, 1.6, 1.8,  1.9,  2.2 /),  &
    GF_OC  = (/ 1.0,  1.2, 1.4, 1.5, 1.6,  1.8,  2.2 /),  &
    GF_EC  = (/ 1.0,  1.0, 1.0, 1.2, 1.4,  1.5,  1.9 /),  &
    GF_SS  = (/ 1.0,  1.6, 1.8, 2.0, 2.4,  2.9,  4.8 /),  &
!.. 550 nm
    Ex_SO4 = (/ 1.114, 1.545, 1.742, 1.862, 2.036, 2.206, 2.558/),  &  !H2SO4
    Ex_OC  = (/ 0.560, 0.652, 0.701, 0.741, 0.821, 0.921, 1.181/),  &
    Ex_EC  = (/ 0.483, 0.484, 0.471, 0.417, 0.363, 0.343, 0.332/),  &
    Ex_SSf = (/ 2.699, 2.547, 2.544, 2.508, 2.444, 2.362, 2.221/),  &
    Ex_SSc = (/ 2.143, 2.103, 2.090, 2.106, 2.084, 2.070, 2.064/)

  logical       :: first_call=.true.
  integer       :: n
  real :: gfSO4, gfOC, gfEC, gfSS, &
    rhoSO4_wet, rhoOC_wet, rhoEC_wet, rhoSS_wet, rhoNO3c_wet,  &
    massGF_SO4, massGF_OC, massGF_EC, massGF_SS,               &
    extSO4, extEC, extOC, extSSf, extSSc, &
    SpecExt_SO4, SpecExt_OC, SpecExt_EC , SpecExt_SSf, SpecExt_SSc, &
    SpecExt_DUf, SpecExt_DUc, SpecExt_NO3f, SpecExt_NO3c, SpecExt_NH4
  real, save :: extDUf=0.0, extDUc=0.0
  real, dimension(NumRH)  :: rh_w
  integer, save :: iDUST_WB_F, iDUST_WB_C

  if(first_call)then
    n=find_index("DUST_WB_F",species(:)%name)
    if(n>0)extDUf=species(n)%ExtC
    iDUST_WB_F=n
    n=find_index("DUST_WB_C",species(:)%name)
    if(n>0)extDUc=species(n)%ExtC
    iDUST_WB_C=n
    first_call=.false.
  endif
  rh_w(:)=0.0
  if(rh<=RelHum(1))then
    rh_w(1)=1.0
  elseif(rh>=RelHum(NumRH))then
    rh_w(NumRH)=1.0
  else
    RHloop: do n = 2, NumRH
      if(rh>RelHum(n)) cycle RHloop
!.. rh interpolation weights
      rh_w(n)  =(rh-RelHum(n-1))/(RelHum(n)-RelHum(n-1))
      rh_w(n-1)=1.0-rh_w(n)
      exit RHloop
    enddo RHloop
    if(debug) &
      write(*,'(a15,i3,5f8.2)') '## Rh >> ', &
        n,rh,RelHum(n),RelHum(n-1),rh_w(n),rh_w(n-1)
  endif      
!.. Interpolate: Growth factors
  gfSO4=dot_product(GF_SO4(:),rh_w(:))
  gfOC =dot_product(GF_OC (:),rh_w(:))
  gfEC =dot_product(GF_EC (:),rh_w(:))
  gfSS =dot_product(GF_SS (:),rh_w(:))
  if(debug)&
    write(*,'(a15,5f8.2)')  '## GFs =  ',rh,gfSO4,gfOC,gfEC,gfSS
!.. Interpolate: Extinction efficiencies
  extSO4=dot_product(Ex_SO4(:),rh_w(:))
  extOC =dot_product(Ex_OC (:),rh_w(:))
  extEC =dot_product(Ex_EC (:),rh_w(:))
  extSSf=dot_product(Ex_SSf(:),rh_w(:))
  extSSc=dot_product(Ex_SSc(:),rh_w(:))
  if(debug)&
    write(*,'(a15,6f10.3)') '## ExtEff ',rh,extSO4,extOC,extEC,extSSf,extSSc 
!.. Density of wet aerosol
!.   rho_w = Vfr_dry*Rho_dry + (1-Vfr_dry)*Rho_water 
!..  where   Vfr_dry = 1/GF**3 (dry volume fraction)
  rhoSO4_wet = rhoSO4 / gfSO4**3 + (1.0-1.0/gfSO4**3) ! *1.0 [g/cm3]
  rhoOC_wet  = rhoOC  / gfOC**3  + (1.0-1.0/gfOC**3)
  rhoEC_wet  = rhoEC  / gfEC**3  + (1.0-1.0/gfEC**3)
  rhoSS_wet  = rhoSS  / gfSS**3  + (1.0-1.0/gfSS**3)
!... Fake
  rhoNO3c_wet= rhoSO4 / gfSS**3  + (1.0-1.0/gfSS**3)
!.. Ratio mass_wet / mass_dry (Mwet/Mdry)
  massGF_SO4 = gfSO4**3 * rhoSO4_wet/rhoSO4
  massGF_OC  = gfOC**3  * rhoOC_wet/rhoOC
  massGF_EC  = gfEC**3  * rhoEC_wet/rhoEC
  massGF_SS  = gfSS**3  * rhoSS_wet/rhoSS
!.. Specific extinction [m2/g] 
!   beta = 3/4 * ExtCoef/rho_wet/rad_eff * Mwet/Mdry
  SpecExt_SO4  = 0.75 * extSO4 * massGF_SO4/ (rhoSO4_wet* Reff_SO4)
  SpecExt_OC   = 0.75 * extOC  * massGF_OC / (rhoOC_wet * Reff_OC)
  SpecExt_EC   = 0.75 * extEC  * massGF_EC / (rhoEC_wet * Reff_EC)
  SpecExt_SSf  = 0.75 * extSSf * massGF_SS / (rhoSS_wet * Reff_SSf)
  SpecExt_SSc  = 0.75 * extSSc * massGF_SS / (rhoSS_wet * Reff_SSc)
if(iDUST_WB_F>0) &
  SpecExt_DUf  = 0.75 * species(iDUST_WB_F)%ExtC/(rhoDU * Reff_DUf)
if(iDUST_WB_C>0) &
  SpecExt_DUc  = 0.75 * species(iDUST_WB_C)%ExtC/(rhoDU * Reff_DUc)
  SpecExt_NH4  = SpecExt_SO4 !... Faking
  SpecExt_NO3f = SpecExt_SO4 !... Faking
!.. VERY CRUDE: Assume NOc sitting on SSc: Q and GF for SSc are applied
  SpecExt_NO3c = 0.75 * extSSc * massGF_SS / (rhoNO3c_wet * Reff_SSc)
!.. Specific extinction [m2/g] 
!   beta = 3/4 * ExtCoef/rho_wet/rad_eff * Mwet/Mdry
!        = 3/4 * ExtCoef/rho_wet/rad_eff * Gf^3*rho_wet/rho_dry
!        = 3/4 * ExtCoef/rad_eff * Gf^3/rho_dry
! SpecExt_SO4  = 0.75 * extSO4/Reff_SO4 * gfSO4**3/rhoSO4
! SpecExt_OC   = 0.75 * extOC /Reff_OC  * gfOC**3 /rhoOC 
! SpecExt_EC   = 0.75 * extEC /Reff_EC  * gfEC**3 /rhoEC 
! SpecExt_SSf  = 0.75 * extSSf/Reff_SSf * gfSS**3 /rhoSS 
! SpecExt_SSc  = 0.75 * extSSc/Reff_SSc * gfSS**3 /rhoSS 
! SpecExt_DUf  = 0.75 * extDUf/Reff_DUf           /rhoDU 
! SpecExt_DUc  = 0.75 * extDUc/Reff_DUc           /rhoDU 
!.. Specific extinction
  extC(:)=species(gtot(:))%ExtC         ! assume dry
  do n = 1, ntot
    select case(species(gtot(n))%name)
      case("SO4"  )   ;extC(n)=SpecExt_SO4
      case("NO3_F")   ;extC(n)=SpecExt_NO3f
      case("NO3_C")   ;extC(n)=SpecExt_NO3f*0.3+SpecExt_NO3c*0.7
      case("NH4_F")   ;extC(n)=SpecExt_NH4
!.. Assume NO3f & NH4f as SO4
!     case("SO4","NO3_F","NH4_F");extC(n)=SpecExt_SO4
!.. Assume NO3c sitting on SSc: Q and GF for SSc are applied
!     case("NO3_C")              ;extC(n)=0.3*SpecExt_SO4 &
!                                        +0.7*SpecExt_SSc*rhoSS/rhoNO3c
      case("EC_F_FFUEL_NEW","EC_F_FFUEL_AGE",&
           "EC_F_WOOD_NEW" ,"EC_F_WOOD_AGE" ,&
           "FFIRE_BC");extC(n)=SpecExt_EC
!     case("EC_C_FFUEL" ,"EC_C_WOOD" );extC(n)=SpecExt_ECc
!     case("POM_F_FFUEL","POM_F_WOOD");extC(n)=SpecExt_OC
!!    case("POM_C_FFUEL"             );extC(n)=SpecExt_OCc
      case("PART_OM_F","FFIRE_OM"    );extC(n)=SpecExt_OC
      case("SEASALT_F");extC(n)=SpecExt_SSf
      case("SEASALT_C");extC(n)=SpecExt_SSc
      case("PPM25","REMPPM25",&
           "PPM25_FIRE","FFIRE_REMPPM25",&
           "DUST_WB_F","DUST_SAH_F");extC(n)=SpecExt_DUf
      case("PPM_C","REMPPM_C",&
           "DUST_WB_C","DUST_SAH_C");extC(n)=SpecExt_DUc
!!    case default     ;extC(n)=0.0     ! disregard un-listed components
    endselect
  enddo
endfunction wetExtC
! <---------------------------------------------------------->
subroutine AOD_Ext (i,j,debug)
!------------------------------------------------
! Calculates AOD
!-------------------------------------------------
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug

  integer, parameter      :: NumAOD=size(AOD_GROUP)
  integer                 :: k, n
  real, dimension(NumAOD) :: kext

!-----------------------------------------------------------------
!   AOD_GROUP = (/ SO4, NO3_F, NH4_F, EC_F_NEW, EC_F_AGE, POC_F, &
!       EXTC  = (/ 8.5, 8.5,   8.5,   7.5,      11.0,     5.7,   &  &.0 for POC!
!                  SEASALT_F, SEASALT_C, DUST_NAT_F, DUST_NAT_C /)
!                  3.0,        0.4       1.0,         0.3,      /)
!__________________________________________________________________

!.. Debug variables for individual components
  logical :: first_call=.true.
  logical, dimension(NumAOD) :: &
    mask_SO4=.false.,mask_NO3=.false.,mask_NH4=.false.,mask_EC=.false.,&
    mask_OM=.false.,mask_SS=.false.,mask_DU=.false.            !mask_POM=.false.

  real :: ext_SO4,ext_NO3,ext_NH4,ext_EC,ext_OM,ext_SS,ext_DU,&!,ext_POM,
          AOD_SO4,AOD_NO3,AOD_NH4,AOD_EC,AOD_OM,AOD_SS,AOD_DU  !,AOD_POM

  if(first_call)then
    do n = 1, NumAOD
      select case(species(AOD_GROUP(n))%name)
        case("SO4")          ;mask_SO4(n)=.true.
        case("NO3_F","NO3_C");mask_NO3(n)=.true.
        case("NH4_F")        ;mask_NH4(n)=.true.
        case("EC_F_FFUEL_NEW","EC_F_FFUEL_AGE",&
             "EC_F_WOOD_NEW" ,"EC_F_WOOD_AGE" ,&
             "FFIRE_BC")     ;mask_EC (n)=.true.
!       case("EC_C_FFUEL","EC_C_WOOD"  );mask_EC (n)=.true.
!       case("POM_F_FFUEL","POM_F_WOOD");mask_POM(n)=.true.
!!      case("POM_C_FFUEL"             );mask_POM(k)=.true.
        case("PART_OM_F","FFIRE_OM");mask_OM(n)=.true.
!       case("POM_C_FFUEL"         );mask_OM(n)=.true.
        case("SEASALT_F","SEASALT_C");mask_SS(n)=.true.
        case("PPM25","PPM25_FIRE","PPM_C",&
             "REMPPM25","FFIRE_REMPPM25","REMPPM_C",&
             "DUST_WB_F","DUST_SAH_F",&
             "DUST_WB_C","DUST_SAH_C");mask_DU(n)=.true.
      endselect
    enddo
    first_call=.false.
  endif

  if(debug)then
    write(*,*) ' #### in AOD module  ###'
    AOD_SO4 = 0.0
    AOD_NO3 = 0.0
    AOD_NH4 = 0.0
    AOD_EC  = 0.0
    AOD_OM  = 0.0
    AOD_SS  = 0.0
    AOD_DU  = 0.0
!   AOD_POM = 0.0
  endif

  AOD(i,j) = 0.0
  kext(:)  = 0.0

  do k = KCHEMTOP, KMAX_MID    !_______________ vertical layer loop
!.. ===========================================================================
!..  Extinction coefficients: Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3]
!..                           summed up for all components  
!..    xn_2d(spec,k)*1.e15 *species(ispec)%molwt/AVOG   [molec/m3] -> [ng/m3]
!..                                                   [ng/m3 ] * 1e-9 -> [g/m3]
!..=>  xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
!.. ===========================================================================

!.. Specific extinction
! Dry extinction coeficients
!   kext(:)=xn_2d(AOD_GROUP,k)*species(AOD_GROUP)%molwt &
!                             *species(AOD_GROUP)%ExtC*1.0e6/AVOG
! Wet extinction coeficients
    kext(:)=xn_2d(AOD_GROUP,k)*species(AOD_GROUP)%molwt*1.0e6/AVOG &
           *wetExtC(NumAOD,AOD_GROUP,rh(k),debug.and.(k==KMAX_MID))

!.. Aerosol extinction optical depth : integral over all vertical layers
!.. [1/m] * [m]
    AOD(i,j) = AOD(i,j) + sum(kext(:)) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))

!.. Extinction coefficient on level k
    Extin_coeff(i,j,k) = sum(kext(:))

    if(debug)then
!.. Extinction coefficients for individual components
      ext_SO4=sum(kext(:),mask_SO4(:))
      ext_NO3=sum(kext(:),mask_NO3(:))
      ext_EC =sum(kext(:),mask_EC (:))
      ext_OM =sum(kext(:),mask_OM (:))
      ext_SS =sum(kext(:),mask_SS (:))
      ext_DU =sum(kext(:),mask_DU (:))
!     ext_POM=sum(kext(:),mask_POM(:))

      if((k==KCHEMTOP+1).or.(k==KMAX_MID))&
        write(*,'(a15,i3,8es10.3)') 'EXTINCs for k =', &
          k, Extin_coeff(i,j,k), &
          ext_SO4,ext_NO3,ext_NH4,ext_EC,ext_OM,ext_SS,ext_DU!,ext_POM

!.. Aerosol optical depth for individual components
      AOD_SO4 = AOD_SO4 + ext_SO4 * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_NO3 = AOD_NO3 + ext_NO3 * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_NH4 = AOD_NH4 + ext_NH4 * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_EC  = AOD_EC  + ext_EC  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_OM  = AOD_OM  + ext_OM  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_SS  = AOD_SS  + ext_SS  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
      AOD_DU  = AOD_DU  + ext_DU  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
!     AOD_POM = AOD_POM + ext_POM * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
    endif
  enddo                        !_______________ vertical layer loop

  if(debug)&
    write(*,'(a30,2i5,8es15.3)') '>>>  AOD / AODs  <<<',   &
      i_fdom(i), j_fdom(j), AOD(i,j),   &
      AOD_SO4,AOD_NO3,AOD_NH4,AOD_EC,AOD_OM,AOD_SS,AOD_DU!,AOD_POM
endsubroutine AOD_Ext
! <---------------------------------------------------------->
! <---------------------------------------------------------->
subroutine AOD_calc (i,j,debug)
!------------------------------------------------
! Calculates AOD.... old, cruder routine
!-------------------------------------------------
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug

  integer :: k
  real :: kext
!-----------------------------------------------------------------
!   AOD_GROUP = (/ SO4, NO3_F, NH4_F, EC_F_NEW, EC_F_AGE, POC_F, &
!       EXTC  = (/ 8.5, 8.5,   8.5,   7.5,      11.0,     5.7,   &
!                  SEASALT_F, SEASALT_C, DUST_NAT_F, DUST_NAT_C /)
!                  3.0,        0.4       1.0,         0.3,      /)
!-----------------------------------------------------------------
  AOD(i,j) = 0.0

  do k = KCHEMTOP, KMAX_MID
!.. ===========================================================================
!..  Extinction coefficients: Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3]
!..                           summed up for all components  
!..    xn_2d(ispec,k)*1.e15 *species(ispec)%molwt/AVOG   [molec/m3] -> [ng/m3]
!..                                                   [ng/m3 ] * 1e-9 -> [g/m3]
!..=>  xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
!.. ===========================================================================
    kext=sum(xn_2d(AOD_GROUP,k)*species(AOD_GROUP)%molwt &
                               *species(AOD_GROUP)%ExtC)*1.0e6/AVOG

!   if(debug.and.(k==18.or.k==KMAX_MID))  &
!     write(*,'(a17,i4,es15.3)')'> Ext. coeff',k,kext

!.. Aerosol extinction optical depth : integral over all vertical layers
!.. [1/m] * [m]
    AOD(i,j) = AOD(i,j) + kext * (z_bnd(i,j,k)-z_bnd(i,j,k+1))

!.. Extinction coefficient on level k
    Extin_coeff(i,j,k) = kext

!   if(debug.and.(k==18.or.k==KMAX_MID))  & 
!     write(*,'(a25,i4,2es15.4,2f8.1)') &
!       '>> Kext AOD for layer',k,kext,AOD(i,j),z_bnd(i,j,k),z_bnd(i,j,k+1)
  enddo

  if(debug) &
    write(*,'(a30,2i4,es15.3)') '>>>  AOD  <<<',i_fdom(i),j_fdom(j),AOD(i,j)
endsubroutine AOD_calc
endmodule AOD_PM_ml

