! < ForestFire_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

module ForestFire_ml
 !----------------------------------------------------------------
 ! Uses emissions from either:
 ! 1) FINNv1 daily data 2002 - 2011
 ! REFERENCE:
 ! Wiedinmyer, C., Akagi, S. K., Yokelson, R. J., Emmons, L. K., Al-Saadi,
 ! J. A., Orlando, J. J., and Soja, A. J.: The Fire INventory from NCAR (FINN) 
 ! - a high resolution global model to estimate the emissions from open 
 !  burning, Geosci. Model Dev. Discuss., 3, 2439-2476, 
 !   doi:10.5194/gmdd-3-2439-2010, 2010.
 ! http://www.geosci-model-dev-discuss.net/3/2439/2010/gmdd-3-2439-2010.html
 !
 ! 2)  GFED 3 (Global Forest Emission database)
 ! http://www.falw.vu/~gwerf/GFED/
 ! Currently programmed for 8-daily data (available for 2001 - 2007)
 !----------------------------------------------------------------
  use CheckStop_ml,      only : CheckStop
  use ChemChemicals_ml,  only : species
  use ChemSpecs_tot_ml
  use EmisDef_ml,        only : ISNAP_NAT &! Fires assigned to SNAP-11 usually
   ,NEMIS_FILE, EMIS_FILE  ! which pollutants are wanted, e.g. sox, pm25

  use GridValues_ml,     only : i_fdom, j_fdom, debug_li, debug_lj, &
                                 debug_proc,xm2,GRIDWIDTH_M
  use Io_ml,             only : PrintLog, datewrite
  use MetFields_ml,      only : z_bnd
  use ModelConstants_ml, only : MasterProc, KMAX_MID, &
                                USE_FOREST_FIRES, DEBUG_FORESTFIRE, &
                                IOU_INST,IOU_HOUR,IOU_HOUR_MEAN, IOU_YEAR
  use NetCDF_ml,         only : ReadField_CDF, Out_netCDF,Real4 ! Reads, writes 
  use OwnDataTypes_ml,   only : Deriv, TXTLEN_SHORT
  use Par_ml,            only : MAXLIMAX, MAXLJMAX, &
                                  me,limax,ljmax
  use PhysicalConstants_ml, only : AVOG
  use Setup_1dfields_ml, only : rcemis
  use SmallUtils_ml,     only : find_index
 ! No. days per year, date-type :
  use TimeDate_ml,only : nydays, nmdays, date, current_date  
implicit none

!  Unimod calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!

  public :: Fire_Emis
  public :: Fire_rcemis
  private :: Export_FireNc

  logical, public, dimension(MAXLIMAX,MAXLJMAX), save ::  burning
  real, private, allocatable, dimension(:,:,:), save :: BiomassBurningEmis

  integer, private, save ::  ieCO  ! index for CO

  character(len=TXTLEN_SHORT), private :: FF_poll
  integer :: iemep



!choose which data set to use (choose only one!)
!GFED  has only 8-daily, FINNv1, but both possibilities are kept for now
    logical :: USE_GFED = .false.
    logical :: USE_FINN = .true.


   !/ Defintions of GFED data. If known, we assign the GFED pollutant which
   !  corresponds to each possible EMEP emission file. Simply add EMEP 
   !  lines as required - be consistent with EmisDefs though. (We can 
   !  have more definitions than used in EmisDefs, but not vice.versa.

   ! Assign mol. wts of the GFED data  where known. If mol. wt set to
   ! zero, the code in Fire-rcemis will use the values from the 
   ! ChemSpecs_ml, species()%molwt.
   ! 
   ! If GFED doesn't have emissionss, set a "-" for GFED, then the
   ! desired emission factor (g/kg DW), and then follow the
   ! example in Fire_setups for NH3:
   !

  ! =======================================================================
   !  Mapping to EMEP species
   ! 
     type, private :: bbtype
       character(len=TXTLEN_SHORT) :: txt
       real :: unitsfac
       real :: frac
       integer :: emep
     end type bbtype

        include 'BiomassBurningMapping.inc'

  ! just to keep track
   real, private, save ::  sum_emis(NBB_DEFS) ! debug
  ! matrix to get from forest-fire species to EMEP ones

    type, private :: bbmap
      integer :: emep
      integer :: bbspec
      real    :: fac
    end type bbmap
    type(bbmap), dimension(NCMSPECS) :: fmap

  ! =======================================================================


contains
  subroutine Fire_Emis(daynumber)

    !.....................................................................
    !**    DESCRIPTION:

    !    Reads forest-fire emissions. So far set up for GFED 8d, but in
    !    principal we can re-code by simply adding alternative
    !    subroutines, e.g. to cope with pre-2001 monthly emissions

    integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


    real    :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    integer :: i,j,nstart, alloc_err
    logical :: my_first_call = .true.   ! DSFF
    logical :: my_first_defs = .true. 
    integer :: dd_old = -1,  n
    real    :: fac, to_kgm2s   
    character(len=TXTLEN_SHORT), dimension(NBBSPECS) :: FF_names

    logical :: calc_remppm = .false.    ! Needed to get REMPPM25
    integer :: ieOC, ieBC, iePM25       !    " / "
    real    :: OMbb
    integer ::  loc_maxemis(2) ! debug

    if( my_first_call ) &
        call PrintLog("Biomass Mapping: "//trim(BiomassBurningMapping),MasterProc)

    nstart = -1 ! reset for GFED
    if(USE_GFED)then
       if (current_date%year<2001) then
          if( my_first_call .and. MasterProc  ) then
             call PrintLog("NO 8d GFED FOREST FIRES BEFORE 2001")
          end if
          my_first_call = .false.
          return
       end if
       if (current_date%year>2007) then
          if( my_first_call .and. MasterProc  ) then
             call PrintLog("NO 8d GFED FOREST FIRES AFTER 2007")
          end if
          my_first_call = .false.
          return
       end if

       if ( DEBUG_FORESTFIRE .and. MasterProc ) then 
          write(*,*) "Into the FIRE days:", current_date%year, &
               daynumber, dd_old, mod ( daynumber, 8 ), my_first_call
       end if

       ! GFED Fire emissions are called at 8 days intervals (1, 9, 17, ....)
       ! 46 values available each year: day 361 is the last one.
       ! Return unless new period

       if ( .not. my_first_call .and. mod ( daynumber, 8 ) /= 1  ) return
       nstart=(current_date%year-2001)*46+(daynumber+7)/8

    endif ! GFED


    if ( DEBUG_FORESTFIRE .and. MasterProc ) then 
       write(*,*) "Into the FIRE days:", current_date%year, &
            daynumber, dd_old, mod ( daynumber, 8 ), my_first_call
    end if
    if (dd_old == daynumber) return   ! Only calculate once per day max
    dd_old = daynumber

    if( my_first_call ) call Fire_setup()   ! Gets InvMolwWtFac

    if(DEBUG_FORESTFIRE .and. MasterProc) &
         write(*,*) "FOREST_FIRE: ", daynumber,nstart, NBBSPECS

    ! We need to look for forest-fire emissions which are equivalent
    ! to the standard emission files:

    ieCO = -999

    do n = 1,  NBBSPECS 

       FF_poll = trim(FF_defs(n)%txt)
       iemep = FF_defs(n)%emep

       if(DEBUG_FORESTFIRE .and. MasterProc) then
          write(*,"( a,2i3,1x,2a8,2i3,a)") "FIRE SETUP: ", &
            n, iemep, trim(species(iemep)%name), &
                trim(FF_poll), len_trim(FF_poll)
       end if

       if(USE_GFED)then
           write(*,*) "WARNING! FFIRE GFED USED! May not be working properly check results!"
             call ReadField_CDF('GLOBAL_ForestFireEmis.nc',FF_poll,&
                  rdemis,nstart,interpol='zero_order',needed=.true.,&
                  UnDef=0.0) ! DS added 20/11 to avoid Ivalues==0 line 2457
             !unit conversion to GFED [g/m2/8day]->[kg/m2/s]
             ! to_kgm2s = 1.0e-3 /(8*24.0*60.0*60.0)
             to_kgm2s = 1.0e-3 /(8*24.0*3600.0)
             rdemis = rdemis * to_kgm2s 
            ! For GFED we have OC and BC as well as PM25. Where the
            ! emep model asks for FFIRE_REMPPM25 we need to get this
            ! by substraction
             if ( trim(species(iemep)%name) == "FFIRE_REMPPM25" ) calc_remppm = .true.
             if ( trim(FF_poll) == "OC" ) ieOC = n
             if ( trim(FF_poll) == "BC" ) ieBC = n
             if ( trim(FF_poll) == "PM25" ) iePM25 = n
      endif
      if(USE_FINN)then
             if( DEBUG_FORESTFIRE ) &
              write(*,*) "FFIRE FINN ", me, n, daynumber,  trim(FF_poll)
             call ReadField_CDF('GLOBAL_ForestFireEmis_FINN.nc',FF_poll,&
                  rdemis,daynumber,interpol='mass_conservative',needed=.true.,&
                  UnDef=0.0, & ! DS added 20/11 to avoid Ivalues==0 line 2457
                  debug_flag=DEBUG_FORESTFIRE)

             !unit conversion now done in BiomassBurningMapping.inc file
             ! since it depends on the chemical scheme
             !unit conversion to FINN 
             !if(FF_poll=='OC'.or.FF_poll=='BC'.or.FF_poll=='PM25')then
             !   !unit is [kg/day]
             !   fac=1.0
             !else
             !   !unit is [mole/day], need kg/day here
             !   fac=species(iemep)%molwt/1000.0
             !endif

              fac = FF_defs(n)%unitsfac

             ![kg or mole /day]->[kg/m2/s]
             fac=fac/(GRIDWIDTH_M*GRIDWIDTH_M * 24*3600)
             do j=1,ljmax
                do i=1,limax
                   rdemis(i,j)=rdemis(i,j)*fac*xm2(i,j)
                enddo
             enddo

            ! For FINN1 we have OC and BC as well as PM25. Where the
            ! emep model asks for FFIRE_REMPPM25 we need to get this
            ! by substraction
             if ( trim(species(iemep)%name) == "FFIRE_REMPPM25" ) calc_remppm = .true. 
             if ( trim(FF_poll) == "OC" ) ieOC = n
             if ( trim(FF_poll) == "BC" ) ieBC = n
             if ( trim(FF_poll) == "PM25" ) iePM25 = n
      endif ! FINN

      if ( my_first_call ) then ! Assume NEMIS_FILE for now
         allocate(BiomassBurningEmis(NBBSPECS,MAXLIMAX,MAXLJMAX),&
              stat=alloc_err)
         call CheckStop( alloc_err, "BB alloc problem")

         my_first_call = .false.

      end if

      !/ CO is special. Keep the index
      if ( trim(FF_poll) == "CO" ) ieCO = n

      ! Assign . units should be [kg/m2/s] here

      BiomassBurningEmis(n,:,:) = rdemis(:,:)

      if ( my_first_defs ) call PrintLog(&
        "ForestFire_ml :: Assigns "//trim(FF_poll) , MasterProc)

      if(DEBUG_FORESTFIRE) sum_emis(n) = sum_emis(n) + &
                               sum( BiomassBurningEmis(n,:,:) )

  end do ! BB_DEFS
      my_first_defs  = .false.

    !/ Logical to let Unimod know if there is any emission here to
    !  worry about

    burning(:,:) =  ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )

    ! If needed, calculate REMPPM25, assuming OM/OC = 1.3  
    !RB: But we assume OM/OC=1.7 for forest fire emissions in the EMEP model? Test changing to 1.7 here also

     if ( calc_remppm ) then
       do j =1, ljmax
       do i =1, limax
        OMbb = 1.7*BiomassBurningEmis(ieOC,i,j) + BiomassBurningEmis(ieBC,i,j) 
        if ( DEBUG_FORESTFIRE .and. MasterProc.and. burning(i,j)  ) then
          write(*,"(a,2i4,3es10.3)") "BURN REMPPM25, ", i_fdom(i), j_fdom(j), &
            BiomassBurningEmis(iePM25,i,j), OMbb, &
            BiomassBurningEmis(iePM25,i,j) - OMbb
        end if
        BiomassBurningEmis(iePM25,i,j) = max( BiomassBurningEmis(iePM25,i,j) -OMbb, 0.0 )
      end do
      end do
     end if

     if ( DEBUG_FORESTFIRE .and. debug_proc ) then
       n = ieCO
       FF_poll = trim(FF_defs(n)%txt)
       loc_maxemis = maxloc(BiomassBurningEmis(n,:,: ) )
       call datewrite("SUM_FF: "//trim(FF_poll),  &
          (/ daynumber, n, i_fdom(loc_maxemis(1)), j_fdom(loc_maxemis(2)) /) ,&
          (/  sum_emis(n), maxval(BiomassBurningEmis(n,:,: ) ), &
             BiomassBurningEmis(n,debug_li,debug_lj) /) &
        ) 
     end if ! debug_proc
  end subroutine Fire_Emis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Fire_setup()

   ! Pre-calculate conversion factors to get from BiomassBurning's kg/m2/s 
   ! to molecules/cm3/s. 
   ! 

    integer :: idef, f, n, alloc_err, iFF, nFF
    integer :: iemep, nemep
    character(len=TXTLEN_SHORT) :: last_FF_poll = "NotSet"
    real :: fac

    nFF   = 0
    nemep = 0

    do n = 1,  NBB_DEFS 

       iemep   = FF_defs(n)%emep
       nemep   = nemep  + 1
       FF_poll = trim(FF_defs(n)%txt)
       if ( FF_poll /= last_FF_poll ) then
          nFF = nFF + 1
          last_FF_poll = trim(FF_poll)
!         if(MasterProc) print *, "MATCH+1 FFIRE ", n, iemep, nemep, nFF, trim(FF_defs(n)%txt)
       end if
      !if(MasterProc) print *, "END MATCH FFIRE ", n, iemep, nemep, nFF, trim(FF_defs(n)%txt)

  !// And one final conversion factor.
  ! The biomassBurning array is kept in kg/m2/s for consistency with other
  ! emissions. We will later convert to molecules/cm3/s after spreading
  ! through a vertical distance dz
  !
  ! If we had E in kg/m2/s, we would then take
  !  E*1.0e3  -> g/m2/s
  !  E*0.1    -> g/cm2/s
  !  E*0.1 /MW * Av -> molec/cm2/s
  !  E*0.001 /MW * Av / DZ -> molec/cm3/s where DZ is spread in m
  !  i.e. fmap should be 0.001*Av/MW
  !  (plus account for the fraction of the inventory assigned to EMEP species)

       fac =  &!BUG  FF_defs(n)%unitsfac  &  ! MW scale if needed
              FF_defs(n)%frac      &  ! mass-fraction split
            * 0.001 * AVOG /species(iemep)%molwt

       fmap(nemep) = bbmap( iemep,  nFF, fac ) 
!       if ( MasterProc ) print *, "FFIRE MAP ", nemep, fmap(nemep)
   end do !  n
  end subroutine Fire_setup
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Fire_rcemis(i,j)

!  Disperses the fire emissions vertically and converts to molecules/cm3/s.

!// Injection height: here over 8 levels. Alternative could be PBL
!   or  equally upto ca. 2*PBL (suggested by Sofiev, GEMS)
!ds QUERY - should the emissions be divided equally by level?
!   - will give a higher mixing ratio for thinner levels

   integer, intent(in) :: i,j

   integer, parameter :: KEMISFIRE = 12
   real, dimension(KEMISFIRE:KMAX_MID) :: invDeltaZfac !  height of layer in m div 9
   integer ::  k, n, iem, iFF

   integer, parameter ::  N_LEVELS = KMAX_MID - KEMISFIRE + 1  ! = 9.0 here


   real    :: origrc, bbe
   logical :: debug_flag


     debug_flag = ( DEBUG_FORESTFIRE .and. &
                     debug_proc .and. i == debug_li .and. j == debug_lj ) 

    if ( debug_flag ) then
      if(BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
        write(*,"(a,5i4,es12.3,f9.3)") "BurningDEBUG ", me, i,j, &
              i_fdom(i), j_fdom(j), BiomassBurningEmis(ieCO,i,j)
     end if

    !/ Here we just divide by the number of levels. Biased towards
    !  different levels since thickness and air content differ. Simple though.

     do k = KEMISFIRE, KMAX_MID
       invDeltaZfac(k) = 1.0/ (z_bnd(i,j,k) - z_bnd(i,j,k+1)) /N_LEVELS
     end do
 
     do n = 1,  NCMSPECS

          iFF = fmap(n)%bbspec
          iem = fmap(n)%emep

          bbe = BiomassBurningEmis(iFF,i,j) * fmap(n)%fac

          origrc = rcemis( fmap(n)%emep, KMAX_MID ) ! just for printout 

           ! distribute vertically:

           do k = KEMISFIRE, KMAX_MID
                rcemis( iem, k ) = rcemis( iem, k )  + bbe * invDeltaZfac(k)
           end do !k

           if ( debug_flag ) then
             k=KMAX_MID
             write(*,"(a,2i3,1x,a8,i4,es10.2,4es10.2)") "FIRERC ",&
              n, IFF, trim(species(fmap(n)%emep)%name), &
                 k, BiomassBurningEmis(iFF,i,j),&
                   fmap(n)%fac, invDeltaZfac(k), origrc, rcemis( iem, k )
           end if

!DSBB    !--  Add up emissions in ktonne ......
!DSBB    !   totemadd(iem) = totemadd(iem) + &
!DSBB    !   tmpemis(iqrc) * dtgrid * xmd(i,j)

        end do ! n
 !       call Export_FireNc() ! Caused problems on last attempt

  end subroutine Fire_rcemis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Export_FireNc()
    type(Deriv) :: def1 ! definition of fields
 
    def1%class='ForestFireEmis' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=1.0      !not used
!FEB2011    def1%inst=.true.      !not used
!FEB2011    def1%year=.false.     !not used
!FEB2011    def1%month=.false.    !not used
!FEB2011    def1%day=.false.      !not used
    def1%name='NOx'        !written
    def1%unit='g/m2'       !written
    def1%name='NOx_zero'       
    def1%name='CO_ASCII'       

    call Out_netCDF(IOU_INST,def1,2,1, BiomassBurningEmis(ieCO,:,:),1.0,&
           CDFtype=Real4,fileName_given='FF.nc')
  end subroutine Export_FireNc

end module ForestFire_ml

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
