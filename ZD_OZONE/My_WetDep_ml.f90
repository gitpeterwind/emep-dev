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
 use ChemSpecs_tot_ml          ! SO2, SO4, etc.
 use ChemSpecs_adv_ml          ! IXADV_SO2, IXADV_SO4, etc.
 use ChemGroups_ml, only : WDEP_OXNGROUP, WDEP_SOXGROUP, WDEP_RDNGROUP
 use ChemChemicals_ml
 use DerivedFields_ml,  only : f_2d, d_2d
 use Io_ml,             only : IO_DEBUG
 use MassBudget_ml,     only : totwdep, wdeploss
 use ModelConstants_ml, only : atwS, atwN, atwPM, IOU_INST, MasterProc, DEBUG_MY_WETDEP
 use OwnDataTypes_ml,   only : depmap
 use SmallUtils_ml,  only : find_index

  implicit none
  private

  public :: Init_WetDep       ! Call from Unimod
  public :: WetDep_Budget     ! Call from Aqueous_ml


  type, public :: WScav
!ss     integer  :: itot     !ds may05 - was adv - confusing
     real  :: W_sca       ! Scavenging ratio/z_Sca/rho = W_sca/1.0e6
     real  :: W_sub       ! same for subcloud
  end type WScav
  

  !dsx integer, public, parameter :: NWETDEP =  14  !SeaS 11  ! Number of solublity classes
!  integer, public, parameter :: NWETDEP =  9  ! Number of solublity classes

  !ds May 2010 New wet-dep system. Pre-define various species or groups
  integer, public, parameter :: NWETDEP_CALC =  8  ! Number of solublity classes

   !  Note - these are for "master" or model species - they do not
   !  need to be present in the chemical scheme. However, the chemical
   !  scheme needs to define wet scavenging after these. If you would
   !  like other characteristics, add them here.

  
  integer, parameter, public :: &
     CWDEP_SO2 = 1,  &
     CWDEP_SO4 = 2,  &
     CWDEP_NH3 = 3,  &
     CWDEP_HNO3 = 4,  &
     CWDEP_H2O2 = 5,  &
     CWDEP_HCHO = 6,  &
     CWDEP_PMf  = 7,  &
     CWDEP_PMc  = 8 

! Import WDepMap array produced from GenChem, with e.g. 
 !integer, public, parameter ::  NWETDEP_ADV  = 14
 ! type(depmap), public, dimension(NWETDEP_ADV), parameter:: WDepMap= (/ &
 !      depmap( IXADV_HNO3, CWDEP_HNO3, -1) &
 !     , depmap( IXADV_HONO, CWDEP_HNO3, -1) &


  include 'CM_WetDep.inc'

 ! And create an array to map from the "calc" to the advected species
 ! Use zeroth column to store number of species in that row

   integer, public, dimension(NWETDEP_CALC,0:NWETDEP_ADV) :: Calc2tot


!EGU - find indices of SO4-like particles (PMc included for now - tmp!)
! for SOA we need a separate array, since Fgas will be needed
!  integer, public, parameter, dimension(5) :: & !EGU:NUM_NONVOLOC+NUM_NONVOLEC) ::&
!   WETDEP_SO4LIKE = (/ aNH4, pNO3_f, PPM25, SSFI, Pb210 /)
!EGU  NONVOLOC, NONVOLEC /)
!    integer, public, parameter, dimension(NUM_NONVOLOC) ::&
!   WETDEP_SO4LIKE = (/ NONVOLOC /)
!  integer, public, parameter, dimension(NUM_VOL) ::&
!   WETDEP_SOALIKE = (/ VOL /)   

  type(WScav), public, dimension(NWETDEP_CALC), save  :: WetDep
  
  integer, public, save  :: WDEP_PREC   ! Used in Aqueous_ml
  integer, private, save :: WDEP_SOX, WDEP_OXN, WDEP_RDN
  integer, private, save :: WDEP_SO2, WDEP_SO4, &
    WDEP_HNO3, WDEP_pNO3_f, WDEP_pNO3_c, WDEP_NH3, WDEP_aNH4

contains

  subroutine Init_WetDep()

    integer :: itot, icalc, n, nc

  !/ INCLOUDFAC is A/v_xmi where A is 5.2 m3 kg-1 s-1, !  and v_xmi is the fallspeed (5 m/s). 
    real, parameter ::  FALLSPEED = 5.0                            ! m/s 
    real, parameter ::  SUBCLFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.1 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.1*SUBCLFAC  & 
                      , EFFCO = 0.4*SUBCLFAC  ! collection efficiency b/clouds - coarse

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * INCLOUDFAC. See also notes in Aqueous_ml.
   !/..                             W_Sca  W_sub
    WetDep(CWDEP_SO2)   = WScav(   0.3,  0.15)   ! Berge+Jakobsen, issh
    WetDep(CWDEP_SO4)   = WScav(   1.0,  EFF25)  ! Berge+Jakobsen, jej
    WetDep(CWDEP_NH3)   = WScav(   1.4,  0.5 )  ! subcloud = 1/3 of cloud for gases
    WetDep(CWDEP_HNO3)  = WScav(   1.4,  0.5)   ! 
    WetDep(CWDEP_H2O2)  = WScav(   1.4,  0.5)   ! 
    WetDep(CWDEP_HCHO)  = WScav(   0.1,  0.03)  ! 
    WetDep(CWDEP_PMf)   = WScav(   1.0,  EFF25) !!
    WetDep(CWDEP_PMc)   = WScav(   1.0,  EFFCO) !!

  ! Other PM compounds treated with SO4LIKE array defined above

   !####################### ds NEW define indices here ##########

     WDEP_PREC= find_index("WDEP_PREC",f_2d(:)%name)
     WDEP_SOX = find_index("WDEP_SOX",f_2d(:)%name)
     WDEP_OXN = find_index("WDEP_OXN",f_2d(:)%name)
     WDEP_RDN = find_index("WDEP_RDN",f_2d(:)%name)

     WDEP_SO2 = find_index("WDEP_SO2",f_2d(:)%name)
     WDEP_SO4 = find_index("WDEP_SO4",f_2d(:)%name)
     WDEP_HNO3 = find_index("WDEP_HNO3",f_2d(:)%name)
     WDEP_pNO3_f = find_index("WDEP_pNO3_f",f_2d(:)%name)
     WDEP_pNO3_c = find_index("WDEP_pNO3_c",f_2d(:)%name)
     WDEP_NH3 = find_index("WDEP_NH3",f_2d(:)%name)
     WDEP_aNH4 = find_index("WDEP_aNH4",f_2d(:)%name)
   !####################### ds END define indices here ##########

   ! Now create table to map calc species to actual advected ones:
     Calc2tot = 0
     do n = 1, NWETDEP_ADV
         icalc = WDepMap(n)%calc
         itot  = WDepMap(n)%ind
         Calc2tot(icalc,0) =  Calc2tot(icalc,0)  + 1
         nc = Calc2tot(icalc,0)
     if( MasterProc ) write(6,"(a,4i5)") "CHECKING WetDep Calc2tot ", n,icalc,itot,nc
         Calc2tot(icalc,nc) = itot
      end do 

     if( MasterProc ) then
      write(*,*) "FINAL WetDep Calc2tot "
      do icalc = 1, NWETDEP_CALC
        write(*,"(i3,i4,15(1x,a))") icalc, Calc2tot(icalc,0 ), &
         ( trim( species(Calc2tot(icalc,nc))%name ), &
             nc= 1, Calc2tot(icalc,0 ))
      end do
     end if

  end subroutine Init_WetDep

  subroutine WetDep_Budget(i,j,invgridarea, debug_flag)
     integer,            intent(in) ::  i,j
     real, intent(in)  :: invgridarea
     logical, intent(in)  :: debug_flag

     real :: wdeps, wdepox, wdepred, wdeppm25, wdeppmco
     real :: fS, fN
     integer :: itot, ispec

     fS = atwS * invgridarea
     fN = atwN * invgridarea

     wdeps = 0.0
     do ispec = 1, size(WDEP_SOXGROUP)
        itot = WDEP_SOXGROUP(ispec)
        wdeps = wdeps + wdeploss(itot)
     end do

     wdepred = 0.0
     do ispec = 1, size(WDEP_RDNGROUP)
        itot = WDEP_RDNGROUP(ispec)
        wdepred = wdepred + wdeploss(itot)
     end do

     wdepox  = 0.0
     do ispec = 1, size(WDEP_OXNGROUP)
        itot = WDEP_OXNGROUP(ispec)
        wdepox  = wdepox  + wdeploss(itot)
     end do

!dsMay2010      wdeps = wdeploss(SO2) + wdeploss(SO4)
!dsMay2010      wdepred = wdeploss(NH3)  + wdeploss(aNH4) 
!dsMay2010      wdepox  = wdeploss(HNO3) + wdeploss(pNO3_f) + wdeploss(pNO3_c)

!dsMay2010      wdeppm25= wdeploss(PPM25) 
!dsMay2010      wdeppmco= wdeploss(PPMco) 

!ds Why mix labels and species like this?
! remove for now
!       totwdep(IXADV_SO4)  = totwdep(IXADV_SO4) + wdeps
!       totwdep(IXADV_HNO3) = totwdep(IXADV_HNO3) + wdepox
!       totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)  + wdepred
!       totwdep(IXADV_PPM25)  = totwdep(IXADV_PPM25)  + wdeppm25
!       totwdep(IXADV_PPMco)  = totwdep(IXADV_PPMco)  + wdeppmco

       d_2d(WDEP_SOX,i,j,IOU_INST) = wdeps * fS 
       d_2d(WDEP_OXN,i,j,IOU_INST) = wdepox * fN 
       d_2d(WDEP_RDN,i,j,IOU_INST) = wdepred * fN 

!dsMay2010 skip for now until
       d_2d(WDEP_SO2,i,j,IOU_INST) = wdeploss(SO2) * fS 
       d_2d(WDEP_SO4,i,j,IOU_INST) = wdeploss(SO4) * fS 
       d_2d(WDEP_NH3,i,j,IOU_INST) = wdeploss(NH3) * fN 
       d_2d(WDEP_aNH4,i,j,IOU_INST) = wdeploss(aNH4) * fN 
       d_2d(WDEP_HNO3,i,j,IOU_INST) = wdeploss(HNO3) * fN 
       d_2d(WDEP_pNO3_f,i,j,IOU_INST) = wdeploss(pNO3_f) * fN 
       d_2d(WDEP_pNO3_c,i,j,IOU_INST) = wdeploss(pNO3_c) * fN 

     if ( DEBUG_MY_WETDEP .and. debug_flag ) then
       write(*,"(a,i4,2es12.4)" )   "DEBUG_WET NH3", WDEP_NH3, wdeploss(NH3), wdeploss(NH3) * fN
       write(*,"(a,i4,2es12.4)" )   "DEBUG_WET ANH4", WDEP_aNH4, wdeploss(aNH4), wdeploss(aNH4) * fN
       write(*,"(a,i4,2es12.4)" )   "DEBUG_WET RDN", WDEP_RDN, wdepred, wdepred* fN
     do ispec = 1, size(WDEP_RDNGROUP)
        itot = WDEP_RDNGROUP(ispec)
        wdepred = wdepred + wdeploss(itot)
         write(*,"(a,2i4,2es12.4)" )   "DEBUG_RDN", ispec, itot, wdeploss(itot)
     end do
     end if

!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdeps     ",i,j,wdeps,wdepred,wdepox!EX,wdeppm25,wdeppmco
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_nh3  ",i,j,d_2d(WDEP_NH3,i,j,IOU_INST), wdeploss(NH3)
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_pNO3_f ",i,j,d_2d(WDEP_pNO3_f,i,j,IOU_INST), wdeploss(pNO3_f)
!write(IO_DEBUG,"(a,2i4,10es12.3)") "wdep_pNO3_c ",i,j,d_2d(WDEP_pNO3_c,i,j,IOU_INST), wdeploss(pNO3_c)


  end subroutine WetDep_Budget
end module My_WetDep_ml
