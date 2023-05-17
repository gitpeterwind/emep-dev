! <BiDir_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.51biDir>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
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
module BiDir_module

  ! Bi-directional  exchange module, following methods developed
  ! by Roy Wichink Kruit (RJWK), Dave Simpson, Sebastiaan Hazelhorst
  ! and DEPAC manual 
  ! Labels : D3.4 means DEPAC chapter 3, eqn 4.
  !          F3   means DEPAC Appendix F, eqn 3.
  !
  ! Refs:
  ! RJWK10: Wichink Kruit, R, van Pul, W, Sauter, F, van den Broek,
  ! M Nemitz, E.  Sutton, M, Krol, M. & Holtslag, A. Modeling the
  ! surface-atmosphere exchange of ammonia Atmos. Environ., 2010, 44,
  ! 945 - 957

  ! RJWK12:  Wichink Kruit, R., J Schaap, M, Sauter, F., J van Zanten,
  ! M. C. & van Pul, W. A. J. Modeling the distribution of ammonia
  ! across Europe including bi-directional surface-atmosphere exchange
  ! Biogeosciences, 2012, 9, 5261-5277

  ! RJWK17: Wichink Kruit, R. J, Aben, J, de Vries, W, Sauter, F, van
  ! der Swaluw, E, van Zanten, M. C. & van Pul, W. A. J., Modelling trends
  ! in ammonia in the Netherlands over the 1990-2014 Atmos. Env., 2017,
  ! 154, 20-30

  ! NOTATION - follows EMEP/DO3SE  conventions
  ! =======================================================
  !  LC = land-cover
  !  R, G = canopy or big-leaf values of resistance/conductance
  !  r, g = leaf-level values of resistance/conductance
  !  ext = external leaf/veg surface, was 'w' in CEH work
  !  gs = ground surface
  !  igs = in-canopy+gs, e.g. Rigs = Rinc + Rgs = DEPAC Rsoil,eff
  !  ns = non-stomatal, EMEP includes gs+ext, ie Gns=Gext+Gigs
  !  w   = water (do not confuse with CEH w!)
  !  Rsur   = 1/(Gsto + Gigs + Gext)  ! 1/R = 1/Rsto + 1/Rns
  
  ! ISSUES
  ! =======================================================
  ! 1. Still need to think more about ground-surface Rgs
  ! 2. And water stuff can be added in 2018. Will try to find
  !   ocean concs NH4

  ! FUTURE; Consider Tleaf calculation. EMEP/LOTOS use Tair so far

  !use Resistances_module, only : RnsA_ACP2012, RextA_DEPAC, RgsA_DEPAC
  implicit none
  private

  ! Main:
  public  :: BiDirXconcs    ! Sets Gammas, X-terms, Rext, for all LC at once
  public  :: BiDirFluxes     ! Depositions and Emissions from specific LC

  public  :: BiDirResistances ! Rext, Gigs,  Rsur
  ! Sea-stuff Jul 2019
  public :: BiDirXwaterEuro     ! From Sebastiaan/Roy work, 2019-
  public :: BiDirXwaterOrig     ! from 2019 work, orig DEPAC

  ! From earlier Resistances_module:
  public :: RnsA_ACP2012  ! Rns(NH3),  EMEP ACP2012 method
  public :: RextA_DEPAC   ! Rext(NH3), Rw in CEH/DEPAC notation
  public :: RgsA_DEPAC    ! DEPAC's Rsoil

  ! for testing
  public  :: self_test       ! Used by mk.testX


  real, private, save :: Xext,Xgs,Xsto  ! hazelhos deleted Xwater
  character(len=20), public :: BiDir_errmsg = 'ok'

 ! We create a new type where we can store a lot of info. Used mainly
 ! by EMEP CTM so far, but aims to be flexible?

 type, public :: BiDir_t

   logical :: EuroXwater = .false.
   logical :: OrigXwater = .false.
   logical :: skipForestDisp   = .false.  !< QUERY. Do we use d for forest?
   character(len=20) :: Method = 'NOTSET'  ! FUTURE
   ! allow for long file names
   character(len=500):: InputFile = 'NOTSET'
   character(len=500) :: InputDir  = 'NOTSET'
  end type BiDir_t
  type(BiDir_t), public, save :: BiDir= BiDir_t()

 type, public :: BiDir2D_t
   real, allocatable, dimension(:,:) :: &
      NH3aLT   &  !long term NH3 (monthly or annual), ug/m3
     ,SO2aLT   &  
     ,NH3inst  & ! NH3 @ 3-4m
     ,NHxEmis  &
     ,aSN      & ! Ratio 
     ,Xtot     & ! Save results
     ,NH3_3m   & ! NH3 @ 3-4m
     ,XH3_3m   & ! NH3 @ 3-4m
     ,NHxDep     !DS MAY2023 If USES%BiDirXwaterOrig:

  !JUL2019 SEA STUFF
  ! If USES%BiDirEuroXwater:
   real, allocatable, dimension(:,:) :: &
      sea_nh4   & ! want umole/L
     ,sea_pH    & !
     ,sea_gridT & !
     ,sea_gridS   !

  end type BiDir2D_t
  type(BiDir2D_t), public, save :: BiDir_2d= BiDir2D_t()

!     tmpNHxDep
! TMP until netcdf read sorted out
!  real, public, allocatable, dimension(:,:), save :: BiDir_NH3acc
!  real, public, allocatable, dimension(:,:), save :: BiDir_SO2acc

contains

 !---------------------------------------------------------------------------
 !---------------------------------------------------------------------------
 !Hazelhos autumn 2019: removed NHxDep, Xwater_in -> Xwater, 
 !Hazelhos 07-02-2020: added water, 
 ! Hazelhos 20-03-2020 removed water.
  subroutine BiDirXconcs(t2C,sstK,Xwater,aSN,NH3aInst,NH3aLT,dbg)

      real, intent(in) :: t2C, sstK ! Temps
      real, intent(inout) :: Xwater ! pre-calculated if monthly data
                                    ! available, otherwise -999
        ! Hazelhos, autumn 2020:  Xwater_in -> Xwater, intent(in) -> intent(inout)
      real, intent(in) :: aSN       !  SO2/NH3
      real, intent(in) :: NH3aInst  ! Ammonia at 4 m , ug/m3, instantaneous
      real, intent(in) :: NH3aLT    ! Ammonia at 4 m , ug/m3, monthly or annual
      logical, intent(in) :: dbg
       !Hazelhos 20-03-2020- logical, intent(in) :: water 
      real :: Gamma_micromet        ! from micro-met at 4m
      real :: Gamma_sto             ! conc. at stomatal interface
      real :: Gamma_ext             ! conc. at leaf-surface-water interface
      real :: Gamma_water           ! seas, rivers, see RJWK12
      real :: Gamma_gs              ! ground-surface
      character(len=*), parameter :: dtxt='BiDirXconc:'
      real :: tmpGamma

      real :: nh3i, Tk, fT, fSST, FaSN  !Hazelhos autumn, 2020: Xwater deleted.

      Tk   = t2C + 273.15

      nh3i      = max(NH3aInst, 1.0e-6 ) ! Prevents numerical problems

      Gamma_ext = 1.84e3*nh3i  * exp( -0.11*t2C ) - 850 !RJWK10 [13],F.1

      tmpGamma = Gamma_ext  !DS for debug print
      FaSN = -999.          !DS for debug print

      if ( Gamma_ext > 0.0 ) then  ! co-deposition, RJWK17 [5], aSN <0.83
        FaSN  = 1.1 - 1.32 * aSN 
        Gamma_ext =  Gamma_ext * max( 0.0, FaSN)
      else
        Gamma_ext = 0.0
      end if
      if ( dbg ) write(*,"(a,f7.2,9es12.3)") dtxt//'BiDirFaSN:',t2C,nh3i, &
               tmpGamma, Gamma_ext, FaSN

      Gamma_gs =  0.0 ! See DEPAC 5 Set zero as too uncertain
      !Gamma_gs =  Gamma_ext ! PRELIM !!!!!!! QUERY !!!!!
                            ! How do we deal with bare-soil etc??
                            ! and surface area?

      Gamma_micromet = 362 * NH3aLT                        !RJWK10 [15a]

      Gamma_sto = 4.7* Gamma_micromet * exp( -0.071*t2C )  !RJWK10 [16],F.3

         !Hazelhos 20-03-2020:
         ! The section below should not be in an if statement. There
         ! should be an Xwater and a Xext, Xigs and Xsto calculated for
         ! every cell.  In the function BiDirFluxes(), the dependency on
         ! LU is used to calculate Xtot and the relative contributions of
         ! each LU type.  Since for LU types that are not water, Xwater is
         ! not used, we do not need to calculate it here.  We changed it,
         ! so that  Xsto, Xgs, Xwater and Xext are calculated for every
         ! grid cell, independent of LU-type.


      !if( Xwater < 0 ) then !Hazelhos 17-12-2019:  This caused the problem
      ! that for land grid cells with values from CMEM, no Xsto and Xgs are
      ! calculated, which are needed further on.
      !Hazelhos 20-03-2020- if( water .EQ. .False. ) then!Hazelhos 07-02-2020:
      ! Fixed the issue described above by having it depend on whether the 
      ! grid cell is water whether Xsto, Xgs and Xwater are calculated.
      
      !Hazelhos 20-03-2020: we outcommented Gamma_water, fSST and Xwater below,
      ! Xwater is already calculated in BiDirXwater(). !DS ONLY if USES%BiDirEuroXwater !
      !Hazelhos 20-03-2020- Gamma_water = 0.0               !TMP for 2017 work

      fT   = 2.75e15/Tk   * exp( -1.04e4/Tk )           ! cf RJWK10 [9]
      !Hazelhos 20-03-2020- fSST = 2.75e15/sstK * exp( -1.04e4/sstK )


      Xext = max(0.0, fT * Gamma_ext)                   !RJWK10 [9],F.2  ug/m3

      if( Xext > 1.0e5 ) then ! TMP checks
          BiDir_errmsg=dtxt//'ERROR? BIGEXT '
          print '(a,f7.2,4es10.2)', BiDir_errmsg, tK, fT, aSN, FaSN, nh3i
          return
          !DS Xext = -1 * Xext ! MAKE MINUS TO SIGNAL PROBLEM
      end if

      Xsto   = fT * Gamma_sto                   !RJWK10 [12=9],F.4 ug/m3

      Xgs    = fT * Gamma_gs                    !PRELIM  ug/m3

      if (dbg) then
        write(*,'(a,5es12.5)') dtxt//'fT, Gamma_ext, Xext,Xw', fT, Gamma_ext, Xext,Xwater
        write(*,'(a,5es12.5)') dtxt// 'Gamma_sto, Xsto, Xgs', Gamma_sto, Xsto, Xgs
      end if


      !Hazelhos 20-03-2020-Xwater =  fSST * Gamma_water  !RJWK12, adjusted, check!  ug/m3
      !Hazelhos 20-03-2020- else
      !  Xwater = Xwater_in  ! Hazelhos, autumn 2020: this was not needed
      ! anymore, since we changed the intent(in) of Xwater to intent(inout)
   
      ! Hazelhos+ 14-02-2020: Two things added here. The first statement was
      ! necessary to cover the grid cells with LU water but without a value in
      ! CMEMS. In the second part, we reset Xsto, Xgs and Xext to 0, as they
      ! can contain values from the previous Land use loop. 

      if ( Xwater < 0 ) Xwater = 0 

     !if ( dbg ) write(*,'(a,3es12.3)') dtxt//' Xsto, Xgs, Xext:  ', Xsto, Xgs, Xext
      !Xsto     = 0 !Hazelhos 20-03-2020: For cells that are partly water,
      ! Xsto Xgs and Xext are set to 0. However, in the LU loop, this means
      ! that Xtot is only dependent on emissions from water in these cells.
      ! Fix by removing this?
     !Xgs      = 0
     !Xext     = 0
   !Hazelhos 20-03-2020- end if ! Cell = Water or Land
   
      if ( dbg ) then
        write(*,'(/,a,2(1x,a,f5.1),4(1x,a,es8.1),3(1x,a,f5.1))') dtxt, &
          'tC:', t2C,' aSN:', aSN, &
          ' Gam: ext', Gamma_ext, 'sto', Gamma_sto, 'gs', Gamma_gs
        write(*,'(a,2(1x,a,f5.1),4(1x,a,es8.1),3(1x,a,f5.1))') dtxt, &
          ' Xwater', Xwater, '  Xs: ext', Xext, 'Xsto', Xsto,'Xgs', Xgs

      end if
  end subroutine BiDirXconcs

 !---------------------------------------------------------------------------
 ! Get ground-surface (eg soil), in-canopy and external-leaf resistances,
 !  along with surface canopy resistance
 ! (NB DEPAC has no Rinc for grass etc ... consider later. For now we stick with
 !  EMEP methods which use Rinc whenever SAI significant)

  subroutine BiDirResistances(SAI,fRH,frozen,water,Rinc,Gsto,&
                                 Gext,Gigs,Rsur,debug)
      real, intent(in) :: SAI         ! surface area index (1-sided, m2/m2)
      real, intent(in) :: fRH         ! RH, fraction
      logical, intent(in) :: frozen, water
      real, intent(in) :: Rinc        ! in.canopy resistance  (s/m)   
      real, intent(in) :: Gsto        ! bulk stom. conductance  (m/s)
      real, intent(out) :: Gext       ! bulk external-leaf conductance (m/s)
      real, intent(out) :: Gigs       ! ground-surface (soil, plus Rinc term) (m/s)
      real, intent(out) :: Rsur       ! canopy resistance  (s/m)   
      logical, intent(in), optional :: debug
      logical :: dbg = .false.
      character(len=*), parameter :: dtxt='BiDirRes:'

      real    :: RgsDry, RgsWet, Rgs, Rext, Gns
 
      if ( present(debug) ) dbg = debug

      if ( dbg ) write(*,'(a, 2es9.2)') dtxt//' Rinc, Gsto,: ', Rinc, Gsto

      call  RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)

  ! QUERY HERE
  ! TMP! DAVE ASSUMPTION. As EMEP doesn't treat dry/wet separately, we assume
  ! increasing wetness after RH 70%. Revise in future; maybe use exponential
  ! (doesn't seem to matter very much as Rext not affected!?)

      Rgs = RgsDry
      if ( fRH > 0.7 ) Rgs = Rgs + ( fRH-0.7) * (RgsWet-RgsDry)
      Gigs = 1/(Rinc + Rgs )

      Rext  = RextA_DEPAC(SAI,fRH,frozen)

      Gext  = 1/Rext
      Gns   = Gext + Gigs
      Rsur = 1/( Gsto + Gns )

      if ( dbg ) write(*,'(a,2f6.2,99g9.2)') dtxt, &
              SAI, fRH, RgsDry, RgsWet, Gigs, Gext, Gns, Rsur
      !DSif ( dbg ) write(*,'(a,2L2)') dtxt//' Water, frozen: ',  water, frozen !Hazelhos added print statement

  end subroutine BiDirResistances

 !---------------------------------------------------------------------------
 ! From EMEP, ACP2012 (which is based upon Nemitz, Smith etc from CEH)

  function RnsA_ACP2012(degC,fRH,aSN) result(Rns)
      real, intent(in) :: degC, fRH   ! Temp (C), RH, fraction
      real, intent(in) :: aSN    ! acidity ratio, molar [SO4]/[NH3]
      real :: F1, F2
      real :: Rns

      F1 = 10.0 * log10(degC+2.0) * exp(100.0*(1.0-fRH)/7.0)
      F2 = 10.0**( 1.6769 - 1.1099 * aSN ) !EMEP's acidity fac

      Rns = min( 1.0/22 * F1 * F2, 200.0 )
      Rns = max( Rns, 10.0 )
  end function RnsA_ACP2012
 !----------------------------------------------------------------------------
 ! Rext as used in DEPAC, from DEPAC Ch.4.
 ! Based upon Sutton & Fowler, 1993, but modified for SAI
 ! DEPAC ch4 notes that this value lower than e.g. Nemitz et al 2001 
 ! due to different pollution climates, and that DEPACE's Xw (Xext) 
 ! accounts for sources in the canopy itself.

  elemental function RextA_DEPAC(SAI,fRH,frozen) result(Rext)
     real, intent(in) :: SAI  ! surface area index, 1-sided, m2/m2
     real, intent(in) :: fRH       ! RH, fraction
     logical, intent(in)  :: frozen
     real, parameter :: SAIHaarweg = 3.5
     real :: Rext
     if ( SAI < 1.0e-6 )  then
       Rext = 99.0e9  ! No uptake since no leaves 
     else if (frozen) then
       Rext = 200.0
     else
       Rext= SAIHaarweg / SAI * 2*exp(100*(1-fRH)/12.0)
     end if
  end function RextA_DEPAC
 !----------------------------------------------------------------------------
 ! 
  subroutine RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)
     !> From DEPAC manual, ch. 5, based upon Erisman et al 1994
     !! Calculate Dry and Wet for later processing. Might make RH-dependent
     !! interpolation rather than simple wet/dry
      logical, intent(in)  :: water
      logical, intent(in)  :: frozen ! to distinguish wet/frozen soil (QUERY ice?)
      real,    intent(out) :: RgsDry, RgsWet
      logical, intent(in), optional  :: dbg

      if ( frozen )  then ! QUERY ICE
        RgsDry = 10000
        RgsWet = 10000
      else
        RgsDry = 1000
        RgsWet = 10
      end if
      if ( water ) RgsDry = RgsWet ! ;-)
  end subroutine RgsA_DEPAC

 !---------------------------------------------------------------------------
  ! BiDirFluxes - calculates effective surface (???) compensation point  Xtot
  ! and associated emission terms.
  ! The original EMEP NH3 resistances used co-dep terms from Nemitz/Smith
  ! For BiDir NH3 we use different Rext terms since Xtot will take
  ! care of the NH3 resistances (rephrase....)
  ! Need to re-calculate Rsur since Rns changing

  subroutine BiDirFluxes(water,SAI,Ve,Gsto,Gext,Gigs,Xwater,Xtot,Xemis,debug)!Hazelhos: added Xwater as an argument

      logical, intent(in) :: water
      real, intent(in) :: SAI         ! surface area index (1-sided, m2/m2)       !Hazelhos: is not used here! Remove?
      real, intent(in) :: Ve          ! Exchange velocity (m/s) =1/(Ra+Rb+Rsur)
      real, intent(in) :: Gsto        ! bulk stom. conductance  (m/s)
      real, intent(in) :: Gext        ! bulk external-leaf conductance (m/s)
      real, intent(in) :: Gigs        ! ground-surface (soil plus Rinc term) (m/s)
      real, intent(in) :: Xwater      ! surface compensation point water (ug/m3)  ! Hazelhos added
      real, intent(out)   :: Xtot     ! effective surface comp. point  (ug/m3)
      real, intent(out)   :: Xemis    !  Emissions from Xtot. We want molec/cm2/s
      logical, intent(in), optional :: debug
      logical :: dbg = .false.
      character(len=*),parameter :: dtxt = 'BiDirFlx:'
      real :: Rsur, Gsur
     ! Units F= ug/m3 * m/s = ug/m2/s.  Want molec/cm2/s
      real, parameter :: AVOG = 6.023e23  ! Avogadros
      real, parameter :: toMoleccm2s = 1.0e-6/17.0*AVOG*1.0e-4

      if ( present(debug) ) dbg = debug

      Gsur   = Gsto + Gigs + Gext
      Rsur   = 1/Gsur 
      if ( dbg ) write(*,'(a,4es9.2)') dtxt//'Gsto, Gigs, Gext, Rsur:', &
              Gsto, Gigs, Gext, Rsur

      if ( water ) then
        Xtot = Xwater
        Xemis = Xtot * Ve * toMoleccm2s   ! emis part of D3.6
        !BUG Xemis = Xwater * Gsur * toMoleccm2s
        if ( dbg ) write(*,'(a,3es9.2,L2)') dtxt//'Water:', Xtot, Xemis, Xwater, water
      else
      
      ! Weighted average surface conc./comp. point; Units: ug/m3
      ! Land-based: (water just had Xtot = Xwater)

       !Xtot=( Rsur/Rext * Xext + Rsur/Rgseff * Xgs + Rsur*Gsto * Xsto ) !D3.4
       ! Use Gext and Gsto to simplify zero conductance cases:

        !Hazelhos 20-03-20: Xext, Xgs and Xsto are apparently 0 for partly
        ! water cells. Results in only emissions from the part that is water. Fix!

        Xtot =  Rsur* ( Gext*Xext + Gigs*Xgs + Gsto*Xsto ) !D3.4

     ! Emissions :-)
     !Hazelhos 10-12-19: QUERY: Xemis is not used. Why calculate? We have
     ! BiDirNHxemissions in output now

        Xemis = Xtot * Ve * toMoleccm2s   ! emis part of D3.6 

        if ( dbg )  then
         write(*,'(a,9(3x,a,es12.3))') dtxt, &
           ' Gext', Gext, ' Xext', Xext, ' Gigs', Gigs, ' Xgs',  Xgs
         write(*,'(a,9(1x,a,es12.3))') dtxt, &
           ' Gsto ', Gsto,' Xsto', Xsto, ' => Xtot:',  Xtot, 'Xemis:', Xemis
        end if ! dbg
      end if 

      if ( Xtot > 1.0e5 ) then ! had some bug somewhere
          BiDir_errmsg=dtxt//'ERROR? BIGXTOT '
          if ( Xtot > 1.0e10 )  BiDir_errmsg=dtxt//'VBIGXTOT:'
          print '(a,f6.2,3(f7.4,es8.1),a,2es9.2)', BiDir_errmsg, SAI, &
              Gext, Xext, Gigs, Xgs, Gsto,Xsto, ' RX:', Rsur, Xtot
          return
!DS           Xtot = -1 * Xtot ! Will trigger warning
!DS           Xemis = -999.0e9

      end if

  end subroutine BiDirFluxes

 !---------------------------------------------------------------------------
 !DSMAY23 - copied code from rv4.15 version
  function BiDirXwaterOrig(i,j,sstK,DepNHx) result (Xwater)
    integer, intent(in) :: i,j
    real, intent(in) :: sstK   ! Temps
    real, intent(in) :: DepNHx    ! DEP-dry+wet, NHx, long-term, kg/ha/yr
    real :: Xwater    ! eq conc air 
    real :: fSST
   ! was based upon DDEP, kg/ha/yr, now have wet too..
   ! as N or NH4?
    real :: Gamma_water  ! conc. at leaf-surface-ater interface


      fSST = 2.75e15/sstK * exp( -1.04e4/sstK )
      Gamma_water = 250.0 * DepNHx  !RJWK12

      Xwater =  fSST * Gamma_water  !RJWK12, adjusted, check!  ug/m3
      if ( Xwater < 0 ) Xwater = 0 
  end function BiDirXwaterOrig
 !---------------------------------------------------------------------------
  function BiDirXwaterEuro(i,j,sstK_nwp,ssNH4,sspH,sstK_monthly,S,water,&
    dbg,method) result(Xwater) !Hazelhos 21-02-2020. Added water boolean and dbg.
   ! cf asman 1994 AE
    integer, intent(in) :: i,j
    real, intent(inout) :: sstK_monthly, ssNH4, sspH, S !Hazelhos 21-02-2020: intent(in) --> intent(inout)
    character(len=*), intent(in) :: method

    logical, intent(in) :: water !Hazelhos 21-02-20: added
    logical, intent(in) :: dbg   !Hazelhos 21-02-20: added
    character(len=*), parameter :: dtxt='BiDirXwater:'

    real, intent(in) :: sstK_nwp ! Temps
!     real, intent(in) :: ssNH4    ! umol/L
!     real, intent(in) :: sspH      !
!     real, intent(in) :: S   ! salinity  promille
     real :: Xwater    ! eq conc air

     real :: IonicStrength   ! ionic strength
     real :: gam_nh4, gam_nh3, Knh4, Hnh3, sstK
     real, parameter :: mwNH3 = 17.0, Rgas = 8.2075e-5 ! atm m3 /mol/K

     if ( method ==  'nwpSST') then
       sstK = sstK_nwp
     else
       sstK = sstK_monthly
     end if

     !if ( dbg ) write(*,'(a,4es12.3,L2)') dtxt//' ssNH4, sspH, sstK, S, water', ssNH4, sspH, sstK, S, water
     !Hazelhos 21-02-2020+: Values are to be based on RWS waterinfo analysis. These are dummy values for tests.
     !Hazelhos 20-03-2020:    We should think of a better general way of parameterizing this. Perhaps by reading a global field?
     !                Also consider splitting the section below up per parameter, not only based on ssNH4.
     if ( water .AND. ( ssNH4 > 1e35 ) ) then !Invalid flag is 1e36. Is there a more elegant way? DS QUERY
       ssNH4 = 0.2
       sspH = 8.2
       !sstK = 288.0! Get temp from grid? --> sstK_nwp is already from the grid. so do nothing  for sstK is fine here.
       S = 0.5
       !if ( dbg ) write(*,'(a,4es12.3,L2)') dtxt//' ssNH4, sspH, sstK, S, water', ssNH4, sspH, sstK, S, water
     end if
     !Hazelhos-

     !Hazelhos 21-02-2020+: Note these functions might not be valid for fresh water. For future work: Reparameterize for fresh water bodies somehow?
     gam_nh4 = 0.883 - 0.0768 * log(S)   !order 0.61 for 35%%

     IonicStrength =  0.00147 + 0.01988 * S + 2.08357e-5 *S*S 

     gam_nh3 = 1 + 0.085 * IonicStrength

     Knh4 = 5.67e-10 * exp( -6286 * (1/sstK - 1/298.15))

     Hnh3 = 56*exp( 4092*(1/sstK - 1/298.15) )  ! M /atm

  !Knh4 * Hnh3 ~ 10 000 from BD

     Xwater =  mwNH3 * ssNH4/ &
           ( Rgas * sstK * Hnh3 * (1/gam_nh3 + ( 10**(-sspH))/ (gam_nh4 * Knh4 )) )
    
  end function BiDirXwaterEuro
 !---------------------------------------------------------------------------
 !
  subroutine self_test()
  ! This self-test is now very long since it merges those from the earlier
  ! Resisistances module with this BiDir_module. Should simplify later,
  ! but self_tests are useful ;-)

   real, parameter :: hVeg=1.0, SAI=4.0, ustar=0.5   ! crop example
   real    :: Ra=log(45.0/(0.1*hVeg))  / ( 0.41*ustar) ! z0=0.1 h, neglect d,L
   real    :: Rb=6.0/ustar, Rinc=14*hVeg/ustar
   real    :: degC, fRH, aSN, sst = 283.0,Gsto=0.02
   real    :: Ve_BD       ! Exchange velocity (m/s) =1/(Ra+Rb+Rsur)
   real    :: nh3i, nh3lt = 10.0   ! ug/m3  inst, long-term NH3 
   real    :: Xtot, Xemis   ! BiDir terms
   real    :: RsurBD, Rsur_emep, Rns_emep, GextBD, GigsBD, Gsur_emep
   real    :: Xwater ! SEA STUFF
   integer :: iRH, iT, inh3
   logical :: water=.false., frozen=.false., dbg=.true.
   real    :: RgsDry, RgsWet, GnsDry, GnsWet


   print *, '========= 1) RESISTANCES - SIMPLE ============================================'

  ! Rgs does not depend on meteo (except wet/dry):
   call  RgsA_DEPAC(water,frozen,RgsDry,RgsWet,dbg)

   aSN = 0.1  ! low SO2/NH3 in eg NL

  ! Print out conductances in cm/s, ie 100/R
   do iRH = 70,100,5
     fRH = 0.01 * iRH
     print *, '====================================================='
     print *, 'RH=', iRH , 'aSN=', aSN
     print '(a,a4,2(2x,4a12))', '> ', 'T', &
        'GnsEMEP', 'GextDEPAC','GnsDEPACdry', 'GnsDEPACwet', &
        'VgEMEP', 'VgDEPACdry', ' VgDEPACwet'

     GextBD= 1/RextA_DEPAC(SAI,fRH,frozen)        ! DEPAC
    ! Ggs = 1/(Rinc + Rgs );  Gns = Gext + Ggs
     GnsDry = GextBD + 1/(Rinc+RgsDry) ! DEPAC
     GnsWet = GextBD + 1/(Rinc+RgsWet) ! DEPAC

     do iT = 0, 30, 10

        degC= real(iT)

        print '(a,i4,2(2x,4f12.1))', '> ', iT,&
           100/RnsA_ACP2012(degC,fRH,aSN), & ! EMEP
           100*GextBD, &                     ! DEPAC, just ext
           100*GnsDry, &                     ! DEPAC dry
           100*GnsWet, &                     ! DEPAC wet
         ! Depositions:
           100/(Ra+Rb+RnsA_ACP2012(degC,fRH,aSN)), & ! EMEP
           100/(Ra+Rb+1/GnsDry), &                   ! DEPAC dry
           100/(Ra+Rb+1/GnsWet)                      ! DEPAC wet
     end do
   end do

   print *, '========= 2) RH, T, NH3i ,tests: ============================================'

   aSN = 0.1
   do iRH = 70,100,5
     fRH = 0.01 * iRH
     write(*,'(9(a6,i5))') "== RH ", iRH , 'h',int(hVeg), &
        'Rinc ', nint(Rinc), ' ==============='

     call BiDirResistances(SAI,fRH,frozen,water,Rinc,Gsto,&
                                 GextBD,GigsBD,RsurBD,debug=.true.)

     Ve_BD = 1/( Ra + Rb + RsurBD )

     do inh3 = 5, 20, 5
       nh3i = real(inh3)
       write(*,'(/,a,9(1x,a,f6.1),/)') '==:', 'NH3inst', nh3i, 'NH3LT', nh3lt,&
          'SAI', SAI, 'h', hVeg, 'Ve(cm/s)', 100*Ve_BD

       do iT = 0, 30, 10
         degC = real(iT)

         ! 3)  Sets Gammas and X values 
         Xwater = 0.1     ! Hazelhos added
         call BiDirXconcs(degC, sst, Xwater, aSN=aSN,  &  ! Hazelhos Xwater_in = 0.1 not permitted anymore
            NH3aInst=nh3i,NH3aLT=nh3lt,  dbg=.true.)      ! Hazelhos: removed DepNHx = 5,Hazelhos 20-03-2020 removed water = .true.

!  subroutine BiDirXconcs(t2C,sstK,Xwater,aSN,NH3aInst,NH3aLT,DepNHx, dbg)

          if ( Xext < 0.0 ) stop 'Xext'

         ! 4)  Gets Xtot and new Rsur :
          call BiDirFluxes(water,SAI,Ve_BD,Gsto,GextBD,GigsBD,Xwater,Xtot,Xemis,dbg)  ! Hazelhos: added Xwater


        ! Compare with EMEP, which is temp dependent
          Rns_emep =  RnsA_ACP2012(degC,fRH,aSN)
          Gsur_emep = 1/Rns_emep + Gsto
          Rsur_emep = 1/Gsur_emep 

          write(*,'(a,2f8.2,5x,a,2f7.2)') 'BDcomp (emep,BD)  Rc:', &
           Rsur_emep, RsurBD,'Vg (cm/s):',  100/(Ra+Rb+Rsur_emep), 100*Ve_BD
  
       end do ! iT
     end do ! NH3
   end do ! RH

   ! Jul 2019 SEA STUFF:
    !call BiDirSea(sstK= 285.0,ssNH4=2.6,sspH=8.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater max', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=2.6,sspH=7.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater max pH', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=0.026,sspH=8.0, S=35.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater low', Xwater
    !call BiDirSea(sstK= 285.0,ssNH4=0.026,sspH=8.0, S=5.0,Xwater=Xwater)
    !write(*,*) 'BDsea Xwater lowS', Xwater


  end subroutine self_test
end module Bidir_module
!-----------------------------------------------------------------------------
! And just to test the above....
!TSTEMX program testr
!TSTEMX   use Bidir_module
!TSTEMX   call self_test()
!TSTEMX end program testr
