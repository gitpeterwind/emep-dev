! <My_Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module specifies the "derived" fields, such as accumulated 
  ! precipitation 
  ! or sulphate, daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as netCDF fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  ! 
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined 
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml. 
  !  
  !   Only text strings used here to define wanted data
  !   All data field characteristics should be defined in Derived_ml, e.g.
  !   in f_2d arrays. 
  !   Derived fields such as d_2d only exist in Derived_ml, so are
  !   accessed here through subroutine calls - using just the (i,j) part
  !   of the bigger d_2d arrays
  !---------------------------------------------------------------------------
 
use CheckStop_ml,  only: CheckStop, StopAll
use GenSpec_adv_ml        ! Use IXADV_ indices...
use GenSpec_shl_ml        ! Use IXSHL_ indices...
use GenSpec_tot_ml,  only : SO2, SO4, HCHO, CH3CHO  &   !  For mol. wts.
                           ,NO2, aNO3, pNO3, HNO3, NH3, aNH4, PM25, PMCO &
                           ,O3, PAN, MPAN, SSfi, SSco  !SS=SeaSalt
use GenChemicals_ml, only : species               !  For mol. wts.
use ModelConstants_ml, only : atwS, atwN, ATWAIR  &
                        , SOURCE_RECEPTOR  &  
                        , KMAX_MID & ! =>  z dimension
                        , PPBINV  &  !   1.0e9
                        , MFAC       ! converts roa (kg/m3 to M, molec/cm3)

use Chemfields_ml, only : xn_adv, xn_shl, cfac
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use Met_ml,        only : z_bnd, roa    ! 6c REM: zeta
use Par_ml,    only: me, MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area 
use SmallUtils_ml,  only : AddArray, LenArray, NOT_SET_STRING, WriteArray, &
                            find_index
use TimeDate_ml,   only : current_date
implicit none
private

 public  :: Init_My_Deriv
 public  :: My_DerivFunc 

 private :: misc_xn   &          ! Miscelleaneous Sums and fractions of xn_adv
           ,pm_calc              ! Miscelleaneous PM's


   character(len=8),  public ,parameter :: model='ZD_OZONE'


   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.
   !  Factor 1.0e6 converts from kg/m2/a to mg/m2/a

  !        We normally distinguish source-receptor (SR) stuff from model
  !        evaluation.  The SR runs should use as few as possible outputs
  !        to keep CPU and disc-requirements down. We define first then the
  !        minimum list of outputs for use in SR, then define an extra list
  !        of parameters needed in model evaluation, or even for the base-case
  !        of SR runs.


  !============ parameters for source-receptor modelling: ===================!

    integer, public, parameter :: MAX_NUM_DERIV2D = 200
    integer, public, parameter :: MAX_NUM_DERIV3D =   5
    integer, public, parameter :: TXTLEN_DERIV =   20
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV2D) :: wanted_deriv2d = NOT_SET_STRING
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV3D) ::  wanted_deriv3d = NOT_SET_STRING

    integer, private, save :: mynum_deriv2d
    integer, private, save :: mynum_deriv3d



 

    character(len=TXTLEN_DERIV), public, parameter, dimension(49) :: &
  D2_SR = (/ &
!
!    Particles: components
       "D2_SO4      ","D2_aNO3     ","D2_pNO3     ","D2_aNH4     " &
      ,"D2_PPM25    ","D2_PPMco    ","D2_PM25_H2O " &
!
!    Particles: sums
      ,"D2_SIA      ","D2_PM25     ","D2_PM10     ","D2_PMco     " &
      ,"D2_SS       ","D2_tNO3     " &
!
!    Ozone and AOTs
      ,"D2_O3       ","D2_MAXO3    " &
      ,"D2_AOT30    ","D2_AOT40    ","D2_AOT60    " &
      ,"D2_AOT30f   ","D2_AOT40f   ","D2_AOT60f   ","D2_AOT40c   " &
      ,"D2_EUAOT30WH","D2_EUAOT30DF","D2_EUAOT40WH","D2_EUAOT40DF" &! NB: do not remove without removing from My_DryDep too
      ,"D2_UNAOT30WH","D2_UNAOT30DF","D2_UNAOT40WH","D2_UNAOT40DF" &! NB: do not remove without removing from My_DryDep too
      ,"D2_MMAOT30WH","D2_MMAOT40WH" &! NB: do not remove without removing from My_DryDep too
      ,"D2_SOMO35   ","D2_SOMO0    " &
!
!    NOy-type sums 
      ,"D2_NO2      ","D2_OXN      ","D2_NOX      ","D2_NOZ      " &
      ,"D2_OX       "  &
!
!    Ecosystem - fluxes: 
 ! NB: do not remove without removing from My_DryDep too 
      ,"D2_AFSTDF0  ","D2_AFSTDF16 ","D2_AFSTBF0  ","D2_AFSTBF16 " &
      ,"D2_AFSTCR0  ","D2_AFSTCR3  ","D2_AFSTCR6  " & !
       ,"D2_O3DF     ","D2_O3WH     " &
!
!    Surface  pressure (for cross section):
      ,"PS          " &
  /)

  !============ Extra parameters for model evaluation: ===================!

    character(len=TXTLEN_DERIV), public, parameter, dimension(16) :: &
  D2_EXTRA = (/ &
       "D2_SO2      ","D2_HNO3     ","D2_NH3      ","D2_VOC      "&
      ,"WDEP_SO2","WDEP_SO4","WDEP_HNO3","WDEP_aNO3", "WDEP_pNO3" & 
      ,"WDEP_NH3", "WDEP_aNH4" & 
      ,"D2_REDN     ","D2_SSfi     ","D2_SSco     ","D2_SNOW","D2_SNratio" &
  /)
!      ,"D2_VddACC", "D2_VddCOA" & ! ECO08


 ! ECO08:
  integer, public, parameter :: &   ! Groups for DDEP and WDEP
    SOX_INDEX = -1, OXN_INDEX = -2, RDN_INDEX = -3
  integer, public, dimension(2) ::  DDEP_SOXGROUP = (/ SO2, SO4 /)
  integer, public, dimension(2) ::  DDEP_RDNGROUP = (/ NH3, aNH4 /)
  integer, public, dimension(6) ::  DDEP_OXNGROUP =  &
                         (/ NO2, HNO3, aNO3, pNO3, PAN, MPAN /)
  integer, public, dimension(6) ::  DDEP_GROUP ! Working array, set
      ! size as max of SOX, OXN, RDN

 ! Ecosystem dep output uses receiver land-cover classes (LCs)
 ! which might include several landuse types, e.g. Conif in 
!  D2_SO2_m2Conif. 


  integer, public, save :: nOutDDep, nOutVg ! ECO08
!hf OutRs, OutRns
  integer, public, save :: nOutRs, nOutRns, nOutGns 


 ! Store indices of ecossytem-species depositions
  type, public :: Dep_type
     character(len=TXTLEN_DERIV) :: name ! e.g. DDEP_SO2_m2Conif
     integer :: LC            ! Index of Receiver land-cover (one 
                              ! of Dep_Receivers)
     integer :: Adv           ! index in xn_adv arrays
     integer :: f2d           ! index in f_2d arrays
     character(len=10) :: txt ! text where needed, e.g. "Conif"
     integer :: atw           ! atomic weight where needed
     character(len=10) :: units ! e.g.  mgN/m2
  end type Dep_type

   !ECO08 - specify some species and land-covers we want to output
   ! depositions for in netcdf files. DDEP_LCS must match one of
   ! the DEP_RECEIVERS  in My_DryDep_ml.
   !
    integer, public, parameter, dimension(13) :: &
      DDEP_SPECS = (/ SOX_INDEX, OXN_INDEX, RDN_INDEX, &
           SO2,  SO4, NH3, aNH4, NO2, HNO3, aNO3, pNO3, PAN, MPAN /)

    character(len=TXTLEN_DERIV), public, parameter, dimension(3) :: &
      DDEP_LCS  = (/ "Grid", "Conif", "Seminat" /)

    integer, public, parameter, dimension(7) :: &
      WDEP_SPECS = (/ SO2,  SO4, aNH4, NH3, aNO3, HNO3, pNO3 /)


  ! Have many combinations: species x ecosystems
  type(Dep_type), public, &
     dimension( size(DDEP_SPECS)*size(DDEP_LCS) ), save :: OutDDep

   !ECO08 - specify some species and land-covers we want to output
   ! dep. velocities for in netcdf files. Set in My_DryDep_ml.

    integer, public, parameter, dimension(6) :: &
      VG_SPECS = (/ O3, NH3, SO2, PM25,  PMCO , HNO3/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(4) :: &
      VG_LCS  = (/ "Grid", "CF", "SNL", "GR" /)
    type(Dep_type), public, &
     dimension( size(VG_SPECS)*size(VG_LCS) ),    save :: OutVg

!hf OutRs
    integer, public, parameter, dimension(2) :: &
      Rs_SPECS = (/ NH3, SO2/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      Rs_LCS  = (/ "Grid", "CF", "SNL", "GR" ,"TC"/)
    type(Dep_type), public, & !surface resistance
     dimension( size(Rs_SPECS)*size(Rs_LCS) ),    save :: OutRs 

    integer, public, parameter, dimension(2) :: &
      Rns_SPECS = (/ NH3, SO2/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      Rns_LCS  = (/ "Grid", "CF", "SNL", "GR" ,"TC"/)
    type(Dep_type), public, & !Non-stomatal resistance
     dimension( size(Rns_SPECS)*size(Rns_LCS) ),    save :: OutRns

    integer, public, parameter, dimension(2) :: &
      Gns_SPECS = (/ NH3, SO2/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      Gns_LCS  = (/ "Grid", "CF", "SNL", "GR" ,"TC"/)
    type(Dep_type), public, & !Non-stomatal conductance
     dimension( size(Gns_SPECS)*size(Gns_LCS) ),    save :: OutGns

!----------------------
! Less often needed:
!exc  "D2_CO     ","D2T_HCHO  ","D2T_CH3CHO","D2_VOC    ",
!exc ,"D2_O3CF   ","D2_O3TC   ","D2_O3GR   ","D2_ACCSU  ",
!"D2_FRNIT  ","D2_MAXOH  ","D2_HMIX   ","D2_HMIX00 ","D2_HMIX12 " &
!exc  "D2_PAN    ","D2_AOT20    " /)

   !======= MY_DERIVED SYSTEM ======================================

  ! use character arrays to specify which outputs are wanted

   character(len=TXTLEN_DERIV), public, parameter, dimension(4) :: &
       WDEP_WANTED = (/ "WDEP_PREC", "WDEP_SOX ", "WDEP_OXN ", &
                      "WDEP_RDN " /)   ! WDEP_PM not used


     character(len=TXTLEN_DERIV), public, parameter, dimension(2) :: &
       D3_WANTED = (/ "D3_O3        ","D3_TH        " /)


    logical, private, parameter :: MY_DEBUG = .false.
    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine Init_My_Deriv()
!hf Rs
    integer :: i, nLC, nVg, nRs, nRns, nGns

   ! Build up the array wanted_deriv2d with the required field names

     call AddArray(WDEP_WANTED, wanted_deriv2d, NOT_SET_STRING)
     call AddArray( D2_SR,  wanted_deriv2d, NOT_SET_STRING)

     if ( .not. SOURCE_RECEPTOR ) then !may want extra?
        call AddArray( D2_EXTRA, wanted_deriv2d, NOT_SET_STRING)
     end if

     !------------- Dry Depositions for d_2d -------------------------
     !ECO08: Add species and ecosystem depositions if wanted:
     ! We find the various combinations of gas-species and ecosystem, 
     ! adding them to the derived-type array OutDDep (e.g. => D2_SO4_m2Conif)

      nLC = 0
      do i = 1, size(DDEP_SPECS)
        do n = 1, size(DDEP_LCS) 

          nLC = nLC + 1

          !CLEAN OutDDep(nLC)%LC  = find_index(DDEP_LCS(n),DDEP_RECEIVERS)
          OutDDep(nLC)%txt  = DDEP_LCS(n)

          if ( DDEP_SPECS(i) > 0 ) then ! normal species, e.g. DDEP_SO2_m2CF
             OutDDep(nLC)%name = "DDEP_"  // &
                 trim( species(DDEP_SPECS(i))%name ) // &
                    "_m2" // trim( DDEP_LCS(n) )

            OutDDep(nLC)%Adv  = DDEP_SPECS(i) - NSPEC_SHL

            if ( species(DDEP_SPECS(i))%sulphurs > 0 ) then
               OutDDep(nLC)%atw  = species(DDEP_SPECS(i))%sulphurs * atwS
               OutDDep(nLC)%units  =  "mgS/m2"

            else if ( species(DDEP_SPECS(i))%nitrogens > 0 ) then
               OutDDep(nLC)%atw  = species(DDEP_SPECS(i))%nitrogens * atwN
               OutDDep(nLC)%units  =  "mgN/m2"
               call CheckStop( species(DDEP_SPECS(i))%sulphurs>0 , &
                " ERROR in DDEP_SPECS: BOTH S AND N!"// &
                   species(DDEP_SPECS(i))%name)

             else
               call StopAll("ERROR: OutDDep atw failure "// &
                   species(DDEP_SPECS(i))%name)
             end if
           else ! GROUP
            OutDDep(nLC)%Adv  = DDEP_SPECS(i) ! e.g. -1 for SOX
            if ( DDEP_SPECS(i) == SOX_INDEX ) then
               OutDDep(nLC)%atw   = atwS
               OutDDep(nLC)%units = "mgS/m2"
               OutDDep(nLC)%name = "DDEP_SOX_m2"//trim(DDEP_LCS(n))
            else if ( DDEP_SPECS(i) == OXN_INDEX ) then
               OutDDep(nLC)%atw   = atwN
               OutDDep(nLC)%units = "mgN/m2"
               OutDDep(nLC)%name = "DDEP_OXN_m2"//trim(DDEP_LCS(n))
            else if ( DDEP_SPECS(i) == RDN_INDEX ) then
               OutDDep(nLC)%atw   = atwN
               OutDDep(nLC)%units = "mgN/m2"
               OutDDep(nLC)%name = "DDEP_RDN_m2"//trim(DDEP_LCS(n))
            end if

           end if


           if(MY_DEBUG .and. me==0)  write(6,"(a,4(a,i4))") &
            "DDEP NEW ", &
             OutDDep(nLC)%name, nLC, "LC ", OutDDep(nLC)%LC, &
            " Adv: ", OutDDep(nLC)%Adv, "ATW ",  OutDDep(nLC)%atw
        end do ! DDEP_SPECS 
      end do ! DDEP_LCS
      nOutDDep = nLC

!nEcoWDep = 0
!do i = 1, size(WDEP_SPECS)
!    EcoWDep(nEcoWDep)%name = "WDEP_"  // trim( species(WDEP_SPECS(i))%name ) )
!end do

    call AddArray( OutDDep(:)%name, wanted_deriv2d, NOT_SET_STRING)

      !------------- Deposition velocities for d_2d -------------------------
      !ECO08: Add species and ecosystem depositions if wanted:
      ! We find the various combinations of gas-species and ecosystem, 
      ! adding them to the derived-type array LCC_DDep (e.g. => VGO3_CF)

      nVg = 0
      do i = 1, size(VG_SPECS)
        do n = 1, size(VG_LCS) 

          nVg = nVg + 1

          OutVg(nVg)%Adv  =   VG_SPECS(i) - NSPEC_SHL
          OutVg(nVg)%txt  =   VG_LCS(n)
          OutVg(nVg)%units  =   "cm/s"
          OutVg(nVg)%name = "VG_"  // trim( species(VG_SPECS(i))%name ) &
               // "_" // trim( VG_LCS(n) )
          !Not yet known: OutVg(nVg)%LC
          if(MY_DEBUG .and. me==0) write(6,*) "VGOUT ", i, n, &
              OutVg(nVg)%name, &
              OutVg(nVg)%Adv,OutVg(nVg)%txt, OutVg(nVg)%units
        end do !n
      end do ! i
      nOutVg = nVg

     call AddArray( OutVg(:)%name, wanted_deriv2d, NOT_SET_STRING)
!hf Rs


      !------------- Surface resistance for d_2d -------------------------
      ! We find the various combinations of gas-species and ecosystem, 
      ! adding them to the derived-type array LCC_DDep (e.g. => RsO3_CF)

      nRs = 0
      do i = 1, size(Rs_SPECS)
        do n = 1, size(Rs_LCS) 

          nRs = nRs + 1

          OutRs(nRs)%Adv  =   Rs_SPECS(i) - NSPEC_SHL
          OutRs(nRs)%txt  =   Rs_LCS(n)
          OutRs(nRs)%units  =   "s/m"
          OutRs(nRs)%name = "Rs_"  // trim( species(Rs_SPECS(i))%name ) &
               // "_" // trim( Rs_LCS(n) )
          !Not yet known: OutRs(nRs)%LC
          if(MY_DEBUG .and. me==0) write(6,*) "RsOUT ", i, n, &
              OutRs(nRs)%name, &
              OutRs(nRs)%Adv,OutRs(nRs)%txt, OutRs(nRs)%units
        end do !n
      end do ! i
      nOutRs = nRs

     call AddArray( OutRs(:)%name, wanted_deriv2d, NOT_SET_STRING)


      !------------- Non stomatal Surface resistance for d_2d ----------
      ! We find the various combinations of gas-species and ecosystem, 
      ! adding them to the derived-type array LCC_DDep (e.g. => RnsO3_CF)

      nRns = 0
      do i = 1, size(Rns_SPECS)
        do n = 1, size(Rns_LCS) 

          nRns = nRns + 1

          OutRns(nRns)%Adv  =   Rns_SPECS(i) - NSPEC_SHL
          OutRns(nRns)%txt  =   Rns_LCS(n)
          OutRns(nRns)%units  =   "s/m"
          OutRns(nRns)%name = "Rns_"  // trim( species(Rns_SPECS(i))%name ) &
               // "_" // trim( Rns_LCS(n) )
          !Not yet known: OutRns(nRns)%LC
          if(MY_DEBUG .and. me==0) write(6,*) "RnsOUT ", i, n, &
              OutRns(nRns)%name, &
              OutRns(nRns)%Adv,OutRns(nRns)%txt, OutRns(nRns)%units
        end do !n
      end do ! i
      nOutRns = nRns

     call AddArray( OutRns(:)%name, wanted_deriv2d, NOT_SET_STRING)



      !------------- Non stomatal Surface conductance for d_2d ----------
      ! We find the various combinations of gas-species and ecosystem, 
      ! adding them to the derived-type array LCC_DDep (e.g. => GnsO3_CF)

      nGns = 0
      do i = 1, size(Gns_SPECS)
        do n = 1, size(Gns_LCS) 

          nGns = nGns + 1

          OutGns(nGns)%Adv  =   Gns_SPECS(i) - NSPEC_SHL
          OutGns(nGns)%txt  =   Gns_LCS(n)
          OutGns(nGns)%units  =   "cm/s"
          OutGns(nGns)%name = "Gns_"  // trim( species(Gns_SPECS(i))%name ) &
               // "_" // trim( Gns_LCS(n) )
          !Not yet known: OutGns(nGns)%LC
          if(MY_DEBUG .and. me==0) write(6,*) "GnsOUT ", i, n, &
              OutGns(nGns)%name, &
              OutGns(nGns)%Adv,OutGns(nGns)%txt, OutGns(nGns)%units
        end do !n
      end do ! i
      nOutGns = nGns

     call AddArray( OutGns(:)%name, wanted_deriv2d, NOT_SET_STRING)




     mynum_deriv2d  = LenArray( wanted_deriv2d, NOT_SET_STRING )

   ! ditto wanted_deriv3d....

     !if ( .not. SOURCE_RECEPTOR ) then
     !   call AddArray( D3_WANTED,  wanted_deriv3d, NOT_SET_STRING)
     !end if
     mynum_deriv3d  = LenArray( wanted_deriv3d, NOT_SET_STRING )


     if ( me == 0 ) then
       write(*,*) "Init_My_Deriv, mynum_deriv2d = ", mynum_deriv2d
       if( MY_DEBUG ) then
           do i = 1, mynum_deriv2d
               write(*,*) "DEBUG DERIV2D ", i, mynum_deriv2d, wanted_deriv2d(i)
           end do
       end if
       call WriteArray(wanted_deriv2d,mynum_deriv2d," Wanted 2d array is")
       write(*,*) "Init_My_Deriv, mynum_deriv3d = ", mynum_deriv3d
       call WriteArray(wanted_deriv3d,mynum_deriv3d," Wanted 3d array is")
     
     end if

  end subroutine Init_My_Deriv
 !=========================================================================
  subroutine My_DerivFunc( e_2d, n, class , timefrac, density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml. For example, we need the 
    ! index for IXADV_O3 for AOTs, and this might not be available in the model
    ! we are running (a PM2.5 model for example), so it is better to define 
    ! this function here.

  real, dimension(:,:), intent(inout) :: e_2d  !  (i,j) 2-d extract of d_2d
  integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
  character(len=*), intent(in)    :: class       ! Class of data
  real, intent(in)    :: timefrac    ! Timestep as frationof hour, dt/3600

  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density     
! density = 1 ( or = roa when unit ug)

  select case ( class )

      case ( "OX", "NOX", "NOZ", "TOXN", "TRDN", "FRNIT", "tNO3 ", "SSalt" )

           call misc_xn( e_2d, n, class, density )

      case ( "SIA", "PM10", "PM25", "PMco" )

          call pm_calc(e_2d, n, class,  density)

      case  default

            print *, "WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
     end select
  


  end subroutine My_DerivFunc
 !=========================================================================

  subroutine pm_calc( pm_2d, n, class, density )

    !/--  calulates PM10 = SIA + PPM2.5 + PPMco

    real, dimension(:,:), intent(inout) :: pm_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  

    select case ( class )

    case ( "SIA" ) 

      forall ( i=1:limax, j=1:ljmax )
        !ds d_2d( n, i,j,IOU_INST) = &
        pm_2d( i,j) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j))& 
         * density(i,j)
      end forall

    case ( "PM25" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt *cfac(IXADV_SSfi,i,j))&  !SeaS
         * density(i,j)
      end forall

    case ( "PMco" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSco,i,j,KMAX_MID) *species(SSco)%molwt *cfac(IXADV_SSco,i,j))&  !SeaS
          * density(i,j)
      end forall

    case ( "PM10" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt*cfac(IXADV_SO4,i,j)   &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt*cfac(IXADV_SSfi,i,j) & !SeaS
         + xn_adv(IXADV_SSco,i,j,KMAX_MID)*species(SSco)%molwt*cfac(IXADV_SSco,i,j))& !SeaS
         * density(i,j)
      end forall

    end select

  end subroutine  pm_calc
 !=========================================================================

!=========================================================================

  subroutine misc_xn( e_2d, n, class, density)
    real, dimension(:,:), intent(inout) :: e_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  
! density = 1 ( or = roa when unit ug)


    !/--  adds up sulphate, nitrate, or whatever is defined

    select case ( class )

    case ( "TOXN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
              * density(i,j)
      end forall


! OX for O3 and NO2 trend studies

    case ( "OX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
                xn_adv(IXADV_O3,i,j,KMAX_MID)  * cfac(IXADV_O3,i,j)   &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) 
      end forall

    case ( "NOX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_NO,i,j,KMAX_MID) &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) &
              ) * density(i,j)
      end forall

    case ( "NOZ" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) &
              + xn_adv(IXADV_PAN,i,j,KMAX_MID) * cfac(IXADV_PAN,i,j) &
              + xn_adv(IXADV_MPAN,i,j,KMAX_MID) * cfac(IXADV_MPAN,i,j) &
              + xn_adv(IXADV_NO3,i,j,KMAX_MID) &
              + 2.0* xn_adv(IXADV_N2O5,i,j,KMAX_MID) &
              + xn_adv(IXADV_ISNI,i,j,KMAX_MID) &
              ) * density(i,j)
      end forall


    case ( "TRDN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
               ( xn_adv(IXADV_NH3,i,j,KMAX_MID) * cfac(IXADV_NH3,i,j)    &
              +  xn_adv(IXADV_aNH4,i,j,KMAX_MID) * cfac(IXADV_aNH4,i,j))  &
               * density(i,j)
      end forall


    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
             ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)  &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
       /max(1E-80, (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j))&
            +  xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)    &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j))
      end forall

    case ( "tNO3" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d(  i,j ) = &
              ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) )&
              * density(i,j)
      end forall

    case ( "SSalt" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_SSfi,i,j,KMAX_MID) * cfac(IXADV_SSfi,i,j) &
              + xn_adv(IXADV_SSco,i,j,KMAX_MID) * cfac(IXADV_SSco,i,j) )&
              * density(i,j)
      end forall

    end select
  end subroutine misc_xn
 !=========================================================================
end module My_Derived_ml
