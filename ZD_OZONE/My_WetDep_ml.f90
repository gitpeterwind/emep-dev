! <My_WetDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module My_WetDep_ml
 use Io_ml,             only : IO_DEBUG
 use MassBudget_ml,     only : totwdep, wdeploss
 use ModelConstants_ml, only : atwS, atwN, atwPM
 use Derived_ml,  only : f_2d, d_2d, IOU_INST
 use SmallUtils_ml,  only : find_index


  use GenSpec_tot_ml          ! SO2, SO4, etc.
  use GenSpec_adv_ml          ! IXADV_SO2, IXADV_SO4, etc.
  implicit none
  private

  public :: Init_WetDep       ! Call from Unimod
  public :: WetDep_Budget     ! Call from Aqueous_ml

  type, public :: WScav
     integer  :: itot     !ds may05 - was adv - confusing
     real  :: W_sca       ! Scavenging ratio/z_Sca/rho = W_sca/1.0e6
     real  :: W_sub       ! same for subcloud
  end type WScav
  

  !dsx integer, public, parameter :: NWETDEP =  14  !SeaS 11  ! Number of solublity classes
  integer, public, parameter :: NWETDEP =  9  ! Number of solublity classes

!EGU - find indices of SO4-like particles (PMc included for now - tmp!)
! for SOA we need a separate array, since Fgas will be needed
  integer, public, parameter, dimension(5) :: & !EGU:NUM_NONVOLOC+NUM_NONVOLEC) ::&
   WETDEP_SO4LIKE = (/ aNH4, aNO3, PPM25, SSFI, Pb210 /)
!EGU  NONVOLOC, NONVOLEC /)
!    integer, public, parameter, dimension(NUM_NONVOLOC) ::&
!   WETDEP_SO4LIKE = (/ NONVOLOC /)
!  integer, public, parameter, dimension(NUM_VOL) ::&
!   WETDEP_SOALIKE = (/ VOL /)   


  type(WScav), public, dimension(NWETDEP), save  :: WetDep
  
  integer, public, save  :: WDEP_PREC   ! Used in Aqueous_ml
  integer, private, save :: WDEP_SOX, WDEP_OXN, WDEP_RDN
  integer, private, save :: WDEP_SO2, WDEP_SO4, &
    WDEP_HNO3, WDEP_aNO3, WDEP_pNO3, WDEP_NH3, WDEP_aNH4

contains

  subroutine Init_WetDep()

  !/ INCLOUDFAC is A/v where A is 5.2 m3 kg-1 s-1, !  and v is the fallspeed (5 m/s). 
    real, parameter ::  FALLSPEED = 5.0                            ! m/s 
    real, parameter ::  SUBCLFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.1 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.1*SUBCLFAC  & 
                      , EFFCO = 0.4*SUBCLFAC  ! collection efficiency b/clouds - coarse

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * INCLOUDFAC. See also notes in Aqueous_ml.

   !/..                        W_Sca  W_sub
  
    WetDep(1)   = WScav(SO2,    0.3,  0.15)   ! Berge+Jakobsen, issh
    WetDep(2)   = WScav(SO4,    1.0,  EFF25)  ! Berge+Jakobsen, jej
    WetDep(3)   = WScav(NH3,    1.4,  0.5 )  ! subcloud = 1/3 of cloud for gases
    WetDep(4)   = WScav(HNO3,   1.4,  0.5)   ! 
    WetDep(5)   = WScav(H2O2,   1.4,  0.5)   ! 
    WetDep(6)   = WScav(HCHO,   0.1,  0.03)  ! 
    WetDep(7)   = WScav(pNO3,   1.0,  EFFCO) !!
    WetDep(8)   = WScav(PPMCO,  1.0,  EFFCO)
    WetDep(9)   = WScav(SSCO,   1.0,  EFFCO)   !SeaS

  ! Other PM compounds treated with SO4LIKE array defined above

   !####################### ds NEW define indices here ##########

     WDEP_PREC= find_index("WDEP_PREC",f_2d(:)%name)
     WDEP_SOX = find_index("WDEP_SOX",f_2d(:)%name)
     WDEP_OXN = find_index("WDEP_OXN",f_2d(:)%name)
     WDEP_RDN = find_index("WDEP_RDN",f_2d(:)%name)

     WDEP_SO2 = find_index("WDEP_SO2",f_2d(:)%name)
     WDEP_SO4 = find_index("WDEP_SO4",f_2d(:)%name)
     WDEP_HNO3 = find_index("WDEP_HNO3",f_2d(:)%name)
     WDEP_aNO3 = find_index("WDEP_aNO3",f_2d(:)%name)
     WDEP_pNO3 = find_index("WDEP_pNO3",f_2d(:)%name)
     WDEP_NH3 = find_index("WDEP_NH3",f_2d(:)%name)
     WDEP_aNH4 = find_index("WDEP_aNH4",f_2d(:)%name)
   !####################### ds END define indices here ##########

  end subroutine Init_WetDep

  subroutine WetDep_Budget(i,j,invgridarea)
     integer,            intent(in) ::  i,j
     real, intent(in)  :: invgridarea
     real :: wdeps, wdepox, wdepred, wdeppm25, wdeppmco
     real :: fS, fN
     fS = atwS * invgridarea
     fN = atwN * invgridarea


      wdeps = wdeploss(SO2) + wdeploss(SO4)

      wdepred = wdeploss(NH3)  + wdeploss(aNH4) 
  
      wdepox  = wdeploss(HNO3) + wdeploss(aNO3) + wdeploss(pNO3)

      wdeppm25= wdeploss(PPM25) 
      wdeppmco= wdeploss(PPMco) 

       totwdep(IXADV_SO4)  = totwdep(IXADV_SO4) + wdeps
       totwdep(IXADV_HNO3) = totwdep(IXADV_HNO3) + wdepox
       totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)  + wdepred
       totwdep(IXADV_PPM25)  = totwdep(IXADV_PPM25)  + wdeppm25
       totwdep(IXADV_PPMco)  = totwdep(IXADV_PPMco)  + wdeppmco

       d_2d(WDEP_SOX,i,j,IOU_INST) = wdeps * fS 
       d_2d(WDEP_OXN,i,j,IOU_INST) = wdepox * fN 
       d_2d(WDEP_RDN,i,j,IOU_INST) = wdepred * fN 

       d_2d(WDEP_SO2,i,j,IOU_INST) = wdeploss(SO2) * fS 
       d_2d(WDEP_SO4,i,j,IOU_INST) = wdeploss(SO4) * fS 
       d_2d(WDEP_NH3,i,j,IOU_INST) = wdeploss(NH3) * fN 
       d_2d(WDEP_aNH4,i,j,IOU_INST) = wdeploss(aNH4) * fN 
       d_2d(WDEP_HNO3,i,j,IOU_INST) = wdeploss(HNO3) * fN 
       d_2d(WDEP_aNO3,i,j,IOU_INST) = wdeploss(aNO3) * fN 
       d_2d(WDEP_pNO3,i,j,IOU_INST) = wdeploss(pNO3) * fN 

!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdeps     ",i,j,wdeps,wdepred,wdepox!EX,wdeppm25,wdeppmco
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_nh3  ",i,j,d_2d(WDEP_NH3,i,j,IOU_INST), wdeploss(NH3)
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_aNO3 ",i,j,d_2d(WDEP_aNO3,i,j,IOU_INST), wdeploss(aNO3)
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_pNO3 ",i,j,d_2d(WDEP_pNO3,i,j,IOU_INST), wdeploss(pNO3)


  end subroutine WetDep_Budget
end module My_WetDep_ml
