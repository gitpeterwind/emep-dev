! <Met_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

module Met_ml

  ! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !_____________________________________________________________________________
!  Subroutines:          Called from:
!    MeteoGridRead        Unimod
!    MeteoRead            Unimod
!    MetModel_LandUse     Unimod
!    metvar(numt)         Unimod
!    metint               PhyChem
!    met_derived          metint    
!    tiphys               metvar
! 
!  metvar - This routines postprocess the meteo fields:
!     ! Unit changes, special definitions etc...
! 
!     e.g. converts wind to u_xmj, v_xmi
! 
!     call tiphys(numt)
! 
!  metint
! 
!     !     this routine does the forward linear stepping of the meteorological
!     !     fields read or derived every 3 hours.
!        q(..1) = q(.., 1) + ...
! 
!     call met_derived - gets u_mid, rho_sruf, ustar_nwp, after each dt_advec
! 
!============================================================================= 

  use BLPhysics_ml,         only : KPBL_MAX &
    ,KPBL_MAX                & ! Max exent of PBL calcs, eg. k=10
                               ! or 1 if no RiB method needed
    ,JericevicRiB_Hmix       & ! TESTING
    ,Venkatram_Hmix          & ! TESTING
    ,Zilitinkevich_Hmix      & ! TESTING
    ,SeibertRiB_Hmix         & ! TESTING
    ,BrostWyngaardKz         & ! TESTING
    ,TI_BLphysics            & ! TESTING or orig
    ,NWP_Kz                  & ! hb 23.02.2010 Kz from meteo 
    ,OBRIAN_Kz               &
    ,SigmaKz_2_m2s

  use CheckStop_ml,         only : CheckStop
  use Functions_ml,         only : Exner_tab, Exner_nd
  use GridValues_ml,        only : xmd, i_fdom, j_fdom, METEOfelt, projection &
       ,gl,gb, gb_glob, gl_glob, MIN_ADVGRIDS   &
       ,Poles, Pole_included, xm_i, xm_j, xm2, sigma_bnd,sigma_mid &
       ,xp, yp, fi, GRIDWIDTH_M,ref_latitude     &
       ,grid_north_pole_latitude,grid_north_pole_longitude &
       ,GlobalPosition,DefGrid,gl_stagg,gb_stagg
  use MetFields_ml    ! Imports many fields, too many for "only"
  use ModelConstants_ml,    only : PASCAL, PT, CLOUDTHRES, METSTEP  &
       ,KMAX_BND,KMAX_MID,NMET &
       ,IIFULLDOM, JJFULLDOM, NPROC  &
       ,MasterProc, DEBUG_MET,DEBUG_i, DEBUG_j, identi, V_RAIN, nmax  &
       ,DEBUG_HMIX, DEBUG_Kz & 
       ,nstep
  use Par_ml           ,    only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX, me  &
       ,limax,ljmax,li0,li1,lj0,lj1  &
       ,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC  &
       ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2  &
       ,MFSIZEINP, IRUNBEG,JRUNBEG, tgi0, tgj0,gi0,gj0  &
       ,MSG_INIT3,MSG_READ4, tlimax, tljmax, parinit
  use PhysicalConstants_ml, only : KARMAN, KAPPA, RGAS_KG, CP, GRAV    &
       ,ROWATER, PI
  use TimeDate_ml,          only : current_date, date,Init_nmdays,nmdays, &
       add_secs,timestamp,&
       make_timestamp, make_current_date, nydays
  use Io_ml ,               only : IO_INFIELD, ios, IO_SNOW, IO_ROUGH, open_file
  use ReadField_ml,         only : ReadField ! reads ascii fields
  use netcdf



  implicit none
  private


  INCLUDE 'mpif.h'
  INTEGER MPISTATUS(MPI_STATUS_SIZE),INFO

  logical, private, save      ::  debug_proc = .false.
  integer, private, save      ::  debug_iloc, debug_jloc  ! local coords


  integer, public :: startdate(4)
  integer, save   :: Nhh &         ! number of field stored per 24 hours
       ,nhour_first&  ! time of the first meteo stored
       ,nrec          ! nrec=record in meteofile, for example
  ! (Nhh=8): 1=00:00 2=03:00 ... 8=21:00
  ! if nhour_first=3 then
  ! 1=03:00 2=06:00...8=24:00

  character (len = 100)        ::  field_not_found='field_not_found'


  public :: MeteoGridRead
  public :: MeteoRead
  public :: MetModel_LandUse
  public :: metvar
  public :: metint
  public :: tiphys

contains



  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine MeteoRead(numt)

    !    the subroutine reads meteorological fields and parameters (every
    !    METSTEP-hours) from NetCDF fields-files, divide the fields into
    !       domains    and sends subfields to the processors


    implicit none

    integer, intent(in):: numt

    character (len = 100), save  ::  meteoname   ! name of the meteofile
    character (len = 100)        ::  namefield & ! name of the requested field
         ,validity    ! field is either instaneous or averaged
    integer ::   ndim,nyear,nmonth,nday,nhour,k
    integer ::   nr   ! Fields are interpolate in
                      ! time (NMET = 2): between nr=1 and nr=2


    type(date)      ::  next_inptime             ! hfTD,addhours_to_input
    type(timestamp) ::  ts_now                   ! time in timestamp format


    real :: nsec                                 ! step in seconds


    nr=2 !set to one only when the first time meteo is read


    if(numt == 1)then !first time meteo is read
       nr = 1
       sdot_at_mid = .false.
       foundustar = .false.
       foundsdot = .false.
       foundSST  = .false.
! hb 23.02.2010 Kz from meteo
       foundKz_met = .false.


       next_inptime = current_date

       ! If origin of meteodomain does not coincide with origin of large domain,
       ! xp and yp should be shifted here, and coordinates must be shifted when
       ! meteofields are read (not yet implemented)



    else

       nsec=METSTEP*3600.0 !from hr to sec
       ts_now = make_timestamp(current_date)
       call add_secs(ts_now,nsec)
       next_inptime=make_current_date(ts_now)


    endif




    nyear=next_inptime%year
    nmonth=next_inptime%month

    nday=next_inptime%day
    nhour=next_inptime%hour

    if(  current_date%month == 1 .and.         &
         current_date%day   == 1 .and.         &
         current_date%hour  == 0 )         &
         call Init_nmdays( current_date )



    if(MasterProc .and. DEBUG_MET) write(6,*) &
         '*** nyear,nmonth,nday,nhour,numt,nmdays2'    &
         ,next_inptime%year,next_inptime%month,next_inptime%day    &
         ,next_inptime%hour,numt,nmdays(2)


    !  Read rec=1 in case 00:00 from 1st January is missing
    if((numt-1)*METSTEP<=nhour_first)nrec=0
    nrec=nrec+1





    if(nrec>Nhh.or.nrec==1) then              ! define a new meteo input file
56     FORMAT(a5,i4.4,i2.2,i2.2,a3)
       write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
       if(MasterProc)write(*,*)'reading ',trim(meteoname)
       nrec = 1
       !could open and close file here instead of in Getmeteofield
    endif


    if(MasterProc .and. DEBUG_MET) write(*,*)'nrec,nhour=',nrec,nhour



    ! 3D fields (i,j,k)
    ndim=3

    !note that u_xmj and v_xmi have dimensions 0:MAXLIJMAX instead of 1:MAXLIJMAX
    !u_xmj(i=0) and v_xmi(j=0) are set in metvar

    namefield='u_wind'
    call Getmeteofield(meteoname,namefield,nrec,ndim,     &
         validity,u_xmj(1:MAXLIMAX,1:MAXLJMAX,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))


    namefield='v_wind'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity,v_xmi(1:MAXLIMAX,1:MAXLJMAX,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='specific_humidity'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, q(:,:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='sigma_dot'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, sdot(:,:,:,nr))
    if(validity==field_not_found)then
       foundsdot = .false.
       if(MasterProc)write(*,*)'WARNING: sigma_dot will be derived from horizontal winds '
    else
       foundsdot = .true.
    endif

    namefield='potential_temperature'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, th(:,:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='precipitation'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, pr(:,:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    pr=max(0.0,pr)  ! AMVB 2009-11-06: positive precipitation


    namefield='3D_cloudcover'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, cc3d(:,:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    if(trim(validity)/='averaged')then
       if(MasterProc)write(*,*)'WARNING: 3D cloud cover is not averaged'
    endif

! hb 23.02.2010 Kz from meteo
    if (NWP_Kz) then
        namefield='eddy_diffusion_coefficient'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             validity, Kz_met(:,:,:,nr))
        if(validity==field_not_found)then
           foundKz_met = .false.
        if(MasterProc)write(*,*)'WARNING: Kz will be derived in model '
        else
           foundKz_met = .true.
        endif
        Kz_met=max(0.0,Kz_met)  ! only positive Kz
    end if
    if( MasterProc .and. DEBUG_Kz)then
          write(6,*)               &
         '*** After Kz', sum(Kz_met(:,:,:,nr)), minval(Kz_met(:,:,:,nr)), maxval(Kz_met(:,:,:,nr)),maxval(Kz_met(:,:,KMAX_BND,nr)), DEBUG_Kz, NWP_Kz, nr, nrec, ndim, namefield

    endif






    ! 2D fields (surface) (i,j)
    ndim=2

    namefield='surface_pressure'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, ps(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='temperature_2m'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, t2_nwp(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

!dsNOV2009
    namefield='relative_humidity_2m'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, rh2m(:,:,nr))
    if(validity==field_not_found)then
        if(MasterProc)write(*,*)'WARNING: relative_humidity_2m not found'
        rh2m(:,:,nr) = -999.9  ! ?
    else
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       rh2m(:,:,nr) = 0.01 * rh2m(:,:,nr)  ! Convert to fraction 
    endif

    namefield='surface_flux_sensible_heat'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, fh(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    if(validity=='averaged')fh(:,:,1)=fh(:,:,nr)

    namefield='surface_flux_latent_heat'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, fl(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    if(validity=='averaged')fl(:,:,1)=fl(:,:,nr)

    namefield='surface_stress'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         validity, tau(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    tau=max(0.0,tau)
    if(validity=='averaged')tau(:,:,1)=tau(:,:,nr)

    namefield='sea_surface_temperature'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        validity, sst(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc)write(*,*)'WARNING: sea_surface_temperature not found '
       foundSST = .false.
    else
       foundSST = .true.
    endif

    namefield='snow_depth'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        validity, sdepth(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc)write(*,*)'WARNING: snow depth not found '
       foundsdepth = .false.
    else
       foundsdepth = .true.
    endif


    namefield='fraction_of_ice' !is really percentage
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        validity, ice(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc)write(*,*)'WARNING: ice coverage (%) not found '
       foundice = .false.
    else
       foundice = .true.
    endif

  end subroutine Meteoread


  subroutine MeteoGridRead(cyclicgrid)

    !   the subroutine reads the grid parameters (projection, resolution etc.)
    !   defined by the meteorological fields
    !

    implicit none

    integer,  intent(out)      :: cyclicgrid
    integer                    :: ndim,nyear,nmonth,nday,nhour,k

    character (len = 100),save :: meteoname !name of the meteofile


    nyear=startdate(1)
    nmonth=startdate(2)
    nday=startdate(3)
    nhour=0
    current_date = date(nyear, nmonth, nday, nhour, 0 )
    call Init_nmdays( current_date )



    !*********initialize grid parameters*********
56  FORMAT(a5,i4.4,i2.2,i2.2,a3)
    write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
    if(MasterProc)write(*,*)'looking for ',trim(meteoname)


    call Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm_i,xm_j,xm2,&
         ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
         ,cyclicgrid)


    if(MasterProc .and. DEBUG_MET)then
       write(*,*)'sigma_mid:',(sigma_mid(k),k=1,20)
       write(*,*)'grid resolution:',GRIDWIDTH_M
       write(*,*)'xcoordinate of North Pole, xp:',xp
       write(*,*)'ycoordinate of North Pole, yp:',yp
       write(*,*)'longitude rotation of grid, fi:',fi
       write(*,*)'true distances latitude, ref_latitude:',ref_latitude
    endif

    call DefGrid()!defines: i_fdom,j_fdom,i_local, j_local,xmd,xm2ji,xmdji,
                  !         sigma_bnd,carea,gbacmax,gbacmin,glacmax,glacmin

  end subroutine Meteogridread


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine metvar(numt)

    ! This routines postprocess the meteo fields:
    ! Unit changes, special definitions etc...


    implicit none

    integer, intent(in):: numt

    !    local

    real,   dimension(KMAX_MID)          ::  prhelp, exf2
    real,   dimension(KMAX_BND)          ::  exf1

    real,   dimension(MAXLJMAX,KMAX_MID) :: usnd   ! send in x
    real,   dimension(MAXLIMAX,KMAX_MID) :: vsnd   ! and in y direction
    real,   dimension(MAXLJMAX,KMAX_MID) :: urcv   ! rcv in x
    real,   dimension(MAXLIMAX,KMAX_MID) :: vrcv   ! and in y direction

    real   bm, cm, dm, divt, x1,x2, xkmin, p1, p2, uvh2, uvhs
    real   ri, z00, a2, cdh, fac, fac2, ro, xkh, dvdz, dvdzm, xlmix
    real   ric, arg, sl2,dz2k,dex12
    real   prhelp_sum,divk(KMAX_MID),sumdiv
    real   inv_METSTEP

    integer :: i, j, k, lx1,lx2, nr,info
    integer request_s,request_n,request_e,request_w
    real ::Ps_extended(0:MAXLIMAX+1,0:MAXLJMAX+1),Pmid,Pu1,Pu2,Pv1,Pv2


    nr = 2
    if (numt == 1) then

       nr = 1

       !-------------------------------------------------------------------
       !  Initialisations:

       call Exner_tab()

       ! Look for processor containing debug coordinates
       debug_iloc    = -999
       debug_jloc    = -999

       do i = 1, limax
          do j = 1, ljmax
!             if (DEBUG_MET .and. &
              if( i_fdom(i) == DEBUG_I .and. j_fdom(j) == DEBUG_J ) then
                debug_proc = .true.
                debug_iloc    = i
                debug_jloc    = j
             end if
          end do
       end do
       if( debug_proc ) write(*,*) "DEBUG EXNER me", me, Exner_nd(99500.0)
       !-------------------------------------------------------------------

    end if


    divt = 1./(3600.0*METSTEP)


    if (neighbor(EAST) .ne. NOPROC) then
       do k = 1,KMAX_MID
          do j = 1,ljmax
             usnd(j,k) = u_xmj(limax,j,k,nr)
          enddo
       enddo
       CALL MPI_ISEND( usnd, 8*MAXLJMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, request_e, INFO)
    endif




    if (neighbor(NORTH) .ne. NOPROC) then
       do k = 1,KMAX_MID
          do i = 1,limax
             vsnd(i,k) = v_xmi(i,ljmax,k,nr)
          enddo
       enddo
       CALL MPI_ISEND( vsnd , 8*MAXLIMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, request_n, INFO)
    endif





    !     receive from WEST neighbor if any

    if (neighbor(WEST) .ne. NOPROC) then

       CALL MPI_RECV( urcv, 8*MAXLJMAX*KMAX_MID, MPI_BYTE, &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
       do k = 1,KMAX_MID
      do j = 1,ljmax
             u_xmj(0,j,k,nr) = urcv(j,k)
          enddo
       enddo

    else

       do k = 1,KMAX_MID
          do j = 1,ljmax
             u_xmj(0,j,k,nr) = u_xmj(1,j,k,nr)
          enddo
       enddo


    endif

    !     receive from SOUTH neighbor if any

    if (neighbor(SOUTH) .ne. NOPROC) then

       CALL MPI_RECV( vrcv, 8*MAXLIMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
       do k = 1,KMAX_MID
          do i = 1,limax
             v_xmi(i,0,k,nr) = vrcv(i,k)
          enddo
       enddo

    else

       if(Poles(2)/=1)  then
          do k = 1,KMAX_MID
             do i = 1,limax
                v_xmi(i,0,k,nr) = v_xmi(i,1,k,nr)
             enddo
          enddo
       else
          !"close" the South pole
          do k = 1,KMAX_MID
             do i = 1,limax
                v_xmi(i,0,k,nr) = 0.0
             enddo
          enddo
       endif

    endif




    if (neighbor(NORTH) == NOPROC.and.Poles(1)==1) then
       !"close" the North pole
       do k = 1,KMAX_MID
          do i = 1,limax
             v_xmi(i,ljmax,k,nr) = 0.0
          enddo
       enddo
    endif

    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_WAIT(request_e, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_n, MPISTATUS, INFO)
    endif


    inv_METSTEP = 1.0/METSTEP

    do j = 1,ljmax
       do i = 1,limax

          !     conversion of pressure from mb to Pascal.

          ps(i,j,nr) = ps(i,j,nr)*PASCAL


          ! surface precipitation, mm/hr

          surface_precip(i,j) = pr(i,j,KMAX_MID) * inv_METSTEP


          rho_surf(i,j)  = ps(i,j,nr)/(RGAS_KG * t2_nwp(i,j,nr) )

          !     For MM5 we get u_xmj*, not tau. Since it seems better to
          !     interpolate tau than u_xmj*  between time-steps we convert

          if ( foundustar) then
             tau(i,j,nr)    = ustar_nwp(i,j)*ustar_nwp(i,j)* rho_surf(i,j)
          end if


          prhelp_sum = 0.0
          prhelp(1) = max(pr(i,j,1),0.)

          prhelp_sum = prhelp_sum + prhelp(1)

          !     pr is 3 hours accumulated precipitation in mm in each
          !     layer summed from above. This is first converted to precipitation
          !     release in each layer.

          do k = 2,KMAX_MID

             prhelp(k) = pr(i,j,k) - pr(i,j,k-1)

          enddo

          !   accumulated deposition over 3 hour interval
          !   k=KMAX_MID now includes accumulated precipitation over all layers
          !   evaporation has been set to zero as it is not accounted for in the
          !   wet deposition


          !   Add up in WetDeposition, to have the prec used in the model

          pr(i,j,:) = prhelp(:)*divt




          !   interpolation of sigma dot for half layers

          if(foundsdot.and.sdot_at_mid)then !pw rv1_9_24
             do k = KMAX_MID,2,-1

                sdot(i,j,k,nr) = sdot(i,j,k-1,nr)            &
                     + (sdot(i,j,k,nr)-sdot(i,j,k-1,nr))   &
                     * (sigma_bnd(k)-sigma_mid(k-1))       &
                     / (sigma_mid(k)-sigma_mid(k-1))

             enddo
          endif

          !    set sdot equal to zero at the top and bottom of atmosphere.

          sdot(i,j,KMAX_BND,nr)=0.0
          sdot(i,j,1,nr)=0.0

          !    conversion from % to fractions (<0,1>) for cloud cover
          !    calculation of cc3dmax (see eulmc.inc) -
          !    maximum of cloud fractions for layers above a given layer

          cc3d(i,j,1) = 0.01 * cc3d(i,j,1)
          cc3dmax(i,j,1) = cc3d(i,j,1)


          lwc(i,j,:)=0.
          do k=2,KMAX_MID
             cc3d(i,j,k) = 0.01 * cc3d(i,j,k)
             cc3dmax(i,j,k) = amax1(cc3dmax(i,j,k-1),cc3d(i,j,k-1))
             lwc(i,j,k)=0.6e-6*cc3d(i,j,k)
          enddo

       enddo
    enddo






    !     lines associated with computation of surface diffusion
    !    coefficient are commented

    xkmin = 1.e-3

    !     derive the meteorological parameters from the basic parameters
    !     read from field-files.

    do j = 1,ljmax
       do i = 1,limax
          p1 = sigma_bnd(KMAX_BND)*(ps(i,j,nr) - PT) + PT

          exf1(KMAX_BND) = CP * Exner_nd(p1)

          z_bnd(i,j,KMAX_BND) = 0.0



          do k = KMAX_MID,1,-1

             !     eddy diffusivity in the surface-layer follows the formulation used
             !     in the nwp-model which is based on Louis (1979), (see mc7e.f).

             !     the shorter loop is the inner loop to save memory. the order
             !     of the do loops will be changed on a vector machine.

             !     exner-function of the half-layers

             p1 = sigma_bnd(k)*(ps(i,j,nr) - PT) + PT


         exf1(k) = CP * Exner_nd( p1 )

             p2 = sigma_mid(k)*(ps(i,j,nr) - PT) + PT

             !     exner-function of the full-levels

         exf2(k) = CP * Exner_nd(p2)

             !     height of the half-layers

             z_bnd(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
                  (exf1(k+1) - exf1(k)))/GRAV


             !     height of the full levels.

             z_mid(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
                  (exf1(k+1) - exf2(k)))/GRAV

             roa(i,j,k,nr) = CP*((ps(i,j,nr) - PT)*sigma_mid(k) + PT)/      &
                  (RGAS_KG*th(i,j,k,nr)*exf2(k))

          enddo  ! k

       enddo
    enddo
    !-----------------------------------------------------------------------

    if( DEBUG_MET .and. debug_proc ) then
     ! Note: example results are given in MetFields
       write(*,*) "DEBUG meIJ" , me, limax, ljmax
       do k = 1, KMAX_MID
          write(6,"(a12,2i3,3f12.4)") "DEBUG_Z",me, k, &
            z_bnd(debug_iloc,debug_jloc,k), z_mid(debug_iloc,debug_jloc,k)!,&
!            zm3d(debug_iloc,debug_jloc,k)
       end do
    end if


    !     Horizontal velocity divided by map-factor.

    do k = 1,KMAX_MID
       do j = 1,ljmax
          do i = 0,limax
             u_xmj(i,j,k,nr) = u_xmj(i,j,k,nr)/xm_j(i,j)

             !divide by the scaling in the perpendicular direction to get effective u_xmj
             !(for conformal projections like Polar Stereo, xm_i and xm_j are equal)

      enddo
       enddo
       do j = 0,ljmax
          do i = 1,limax
             v_xmi(i,j,k,nr) = v_xmi(i,j,k,nr)/xm_i(i,j)

             !  divide by the scaling in the perpendicular direction to get effective v_xmi
          enddo
       enddo
    enddo



    if(.not.foundsdot)then
       ! sdot derived from divergence=0 principle

       do j = 1,ljmax
          do i = 1,limax
             Ps_extended(i,j) = Ps(i,j,nr)
          enddo
       enddo
!Get Ps at edges from neighbors
!we reuse usnd, vsnd etc
    if (neighbor(EAST) .ne. NOPROC) then
          do j = 1,ljmax
             usnd(j,1) = ps(limax,j,nr)
          enddo
       CALL MPI_ISEND( usnd, 8*MAXLJMAX, MPI_BYTE,  &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, request_e, INFO)
    endif
    if (neighbor(NORTH) .ne. NOPROC) then
          do i = 1,limax
             vsnd(i,1) = ps(i,ljmax,nr)
          enddo
       CALL MPI_ISEND( vsnd , 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, request_n, INFO)
    endif

    !     receive from WEST neighbor if any
    if (neighbor(WEST) .ne. NOPROC) then
       CALL MPI_RECV( urcv, 8*MAXLJMAX, MPI_BYTE, &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
      do j = 1,ljmax
             Ps_extended(0,j) = urcv(j,1)
          enddo
    else
          do j = 1,ljmax
             Ps_extended(0,j) = Ps_extended(1,j)
          enddo
    endif
    !     receive from SOUTH neighbor if any
    if (neighbor(SOUTH) .ne. NOPROC) then
       CALL MPI_RECV( vrcv, 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
          do i = 1,limax
             Ps_extended(i,0) = vrcv(i,1)
          enddo
    else
          do i = 1,limax
             Ps_extended(i,0) = Ps_extended(i,1)
          enddo
    endif
    if (neighbor(WEST) .ne. NOPROC) then
          do j = 1,ljmax
             usnd(j,2) = ps(1,j,nr)
          enddo
       CALL MPI_ISEND( usnd(1,2), 8*MAXLJMAX, MPI_BYTE,  &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, request_w, INFO)
    endif
    if (neighbor(SOUTH) .ne. NOPROC) then
          do i = 1,limax
             vsnd(i,2) = ps(i,1,nr)
          enddo
       CALL MPI_ISEND( vsnd(1,2) , 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, request_s, INFO)
    endif



    !     receive from EAST neighbor if any
    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_RECV( urcv, 8*MAXLJMAX, MPI_BYTE, &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
      do j = 1,ljmax
             Ps_extended(limax+1,j) = urcv(j,1)
          enddo
    else
          do j = 1,ljmax
             Ps_extended(limax+1,j) = Ps_extended(limax,j)
          enddo
    endif
    !     receive from NORTH neighbor if any
    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_RECV( vrcv, 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
          do i = 1,limax
             Ps_extended(i,ljmax+1) = vrcv(i,1)
          enddo
    else
          do i = 1,limax
             Ps_extended(i,ljmax+1) = Ps_extended(i,ljmax)
          enddo
    endif

    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_WAIT(request_e, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_n, MPISTATUS, INFO)
    endif
    if (neighbor(WEST) .ne. NOPROC) then
       CALL MPI_WAIT(request_w, MPISTATUS, INFO)
    endif

    if (neighbor(SOUTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_s, MPISTATUS, INFO)
    endif


!note that u_xmj and v_xmi have already been divided by xm here
       do j = 1,ljmax
          do i = 1,limax
             Pmid=Ps_extended(i,j)-PT
             Pu1=0.5*(Ps_extended(i-1,j)+Ps_extended(i,j))-PT
             Pu2=0.5*(Ps_extended(i+1,j)+Ps_extended(i,j))-PT
             Pv1=0.5*(Ps_extended(i,j-1)+Ps_extended(i,j))-PT
             Pv2=0.5*(Ps_extended(i,j+1)+Ps_extended(i,j))-PT

             sdot(i,j,KMAX_BND,nr)=0.0
             sdot(i,j,1,nr)=0.0
             sumdiv=0.0
             do k=1,KMAX_MID
                divk(k)=((u_xmj(i,j,k,nr)*Pu2-u_xmj(i-1,j,k,nr)*Pu1)            &
                     + (v_xmi(i,j,k,nr)*Pv2-v_xmi(i,j-1,k,nr)*Pv1))            &
                     * xm2(i,j)*(sigma_bnd(k+1)-sigma_bnd(k))  &
                     / GRIDWIDTH_M/Pmid
                sumdiv=sumdiv+divk(k)
             enddo
             sdot(i,j,KMAX_MID,nr)=-(sigma_bnd(KMAX_MID+1)-sigma_bnd(KMAX_MID))*sumdiv+divk(KMAX_MID)
             do k=KMAX_MID-1,2,-1
                sdot(i,j,k,nr)=sdot(i,j,k+1,nr)-(sigma_bnd(k+1)-sigma_bnd(k))*sumdiv+divk(k)
             enddo
          enddo
       enddo
    endif



    call met_derived !compute derived meteo fields

    call tiphys(numt)

  end subroutine metvar

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






  subroutine metint

    !     this routine does the forward linear stepping of the meteorological
    !     fields read or derived every 3 hours.


    implicit none

    integer :: i,j
    real :: div,ii

    if (nstep.lt.nmax) then

       div = 1./real(nmax-(nstep-1))

       u_xmj(:,:,:,1)    = u_xmj(:,:,:,1)         		&
            + (u_xmj(:,:,:,2) - u_xmj(:,:,:,1))*div
       v_xmi(:,:,:,1)    = v_xmi(:,:,:,1)         		&
            + (v_xmi(:,:,:,2) - v_xmi(:,:,:,1))*div
       sdot(:,:,:,1) = sdot(:,:,:,1)         	&
            + (sdot(:,:,:,2) - sdot(:,:,:,1))*div
       th(:,:,:,1)   = th(:,:,:,1)         		&
            + (th(:,:,:,2) - th(:,:,:,1))*div
       q(:,:,:,1)    = q(:,:,:,1)         		&
            + (q(:,:,:,2) - q(:,:,:,1))*div
       !      ccc(:,:,:,1) = ccc(:,:,:,1)                           & !ASSYCON
       !                   + (ccc(:,:,:,2) - ccc(:,:,:,1))*div       !ASSYCON
       SigmaKz(:,:,:,1)  = SigmaKz(:,:,:,1)         		&
            + (SigmaKz(:,:,:,2) - SigmaKz(:,:,:,1))*div
       roa(:,:,:,1)  = roa(:,:,:,1)         		&
            + (roa(:,:,:,2) - roa(:,:,:,1))*div
       ps(:,:,1)     = ps(:,:,1)         		&
            + (ps(:,:,2) - ps(:,:,1))*div
       t2_nwp(:,:,1) = t2_nwp(:,:,1)         		&
            + (t2_nwp(:,:,2) - t2_nwp(:,:,1))*div
       rh2m(:,:,1) = rh2m(:,:,1)  &
            + (rh2m(:,:,2) - rh2m(:,:,1))*div
       SoilWater(:,:,1) = SoilWater(:,:,1)   &
            + (SoilWater(:,:,2) - SoilWater(:,:,1))*div
       SoilWater_deep(:,:,1) = SoilWater_deep(:,:,1)    &
            + (SoilWater_deep(:,:,2) - SoilWater_deep(:,:,1))*div


       fh(:,:,1)     = fh(:,:,1)         		&
            + (fh(:,:,2) - fh(:,:,1))*div
       fl(:,:,1)     = fl(:,:,1)         		&
            + (fl(:,:,2) - fl(:,:,1))*div
       tau(:,:,1)    = tau(:,:,1)         		&
            + (tau(:,:,2) - tau(:,:,1))*div
       sst(:,:,1)    = sst(:,:,1)         		&
            + (sst(:,:,2)   - sst(:,:,1))*div
       sdepth(:,:,1)    = sdepth(:,:,1)         		&
            + (sdepth(:,:,2)   - sdepth(:,:,1))*div
       ice(:,:,1)    = ice(:,:,1)         		&
            + (ice(:,:,2)   - ice(:,:,1))*div

       !  precipitation and cloud cover are no longer interpolated

    else

       !     assign the the meteorological data at time-level 2 to level 1 for
       !     the next 6 hours integration period before leaving the inner loop.

       u_xmj(:,:,:,1)    = u_xmj(:,:,:,2)
       v_xmi(:,:,:,1)    = v_xmi(:,:,:,2)
       sdot(:,:,:,1) = sdot(:,:,:,2)
       th(:,:,:,1)   = th(:,:,:,2)
       q(:,:,:,1)    = q(:,:,:,2)
       !       ccc(:,:,:,1) = ccc(:,:,:,2)   !ASSYCON
       SigmaKz(:,:,:,1)  = SigmaKz(:,:,:,2)
       roa(:,:,:,1)  = roa(:,:,:,2)
       !  - note we need pressure first before surface_pressure
       ps(:,:,1)     = ps(:,:,2)
       t2_nwp(:,:,1) = t2_nwp(:,:,2)
       rh2m(:,:,1) = rh2m(:,:,2)
       SoilWater(:,:,1) = SoilWater(:,:,2)
       SoilWater_deep(:,:,1) = SoilWater_deep(:,:,2)
       sdepth(:,:,1) = sdepth(:,:,2)
       ice(:,:,1) = ice(:,:,2)

       fh(:,:,1)     = fh(:,:,2)
       tau(:,:,1)    = tau(:,:,2)
       fl(:,:,1)     = fl(:,:,2)

       sst(:,:,1)    = sst(:,:,2)

    endif

    call met_derived !update derived meteo fields

  end subroutine metint

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine met_derived

    ! This routine calculates fields derived from meteofields.
    ! The interpolation in time is done for the meteofields and the
    ! fields here are derived from the interpolated fields after
    ! each interpolation (i.e. every dt_advec).
    ! CPU costly fields (those with special functions like log )
    ! can be computed in metvar only once every METSTEP and interpolated
    ! in metint.

    !horizontal wind speed (averaged over the four edges)
    !Note that u_xmj and v_xmi are wind velocities divided by xm
    !At present u_ref is defined at KMAX_MID


    implicit none
    integer ::i,j, k
    logical :: DEBUG_DERIV = .false.

    do k = KPBL_MAX, KMAX_MID
    do j = 1,ljmax
       do i = 1,limax
           u_mid(i,j,k) = 0.5*( u_xmj(i,j,k,1)*xm_j(i,j) + &
                                u_xmj(i-1,j,k,1)*xm_j(i-1,j) )
           v_mid(i,j,k) = 0.5*( v_xmi(i,j-1,k,1)*xm_i(i,j-1) + &
                                v_xmi(i,j,k,1)*xm_i(i,j))
       enddo
    enddo
    enddo !k
    do j = 1,ljmax
       do i = 1,limax
          u_ref(i,j)= sqrt( u_mid(i,j,KMAX_MID)**2 + v_mid(i,j,KMAX_MID)**2 )
       enddo
    enddo

    ! Tmp ustar solution. May need re-consideration for MM5 etc., but
    ! basic principal should be that fm is interpolated with time, and
    ! ustar derived from this.

    !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

    forall( i=1:limax, j=1:ljmax )
       rho_surf(i,j)  = ps(i,j,1)/(RGAS_KG * t2_nwp(i,j,1) )
    end forall

    if(.not. foundustar)then
       forall( i=1:limax, j=1:ljmax )
          ustar_nwp(i,j)   = sqrt( tau(i,j,1)/rho_surf(i,j) )
       end forall
    endif


    forall( i=1:limax, j=1:ljmax )
       ustar_nwp(i,j) = max( ustar_nwp(i,j), 1.0e-5 )
    end forall

    if ( DEBUG_DERIV .and. debug_proc ) then
       i = debug_iloc
       j = debug_jloc
       write(*,*) "MET_DERIV DONE ", me, ustar_nwp(i,j), rho_surf(i,j)
    end if
    !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


  end subroutine met_derived

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




  subroutine MetModel_LandUse(callnum)

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !     This subroutine reads parameterfields from file
    !     reading surface roughness classes from file: rough.dat
    !     reading snow                      from file: snowc.dat
    !
    !     ... fields as used in meteorological model
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    implicit none

    integer, intent(in) :: callnum
    integer ::  i,j, err

    real, dimension(MAXLIMAX,MAXLJMAX) :: r_class  ! Roughness (real)

    character*20 fname
    logical :: needed_found

    ios = 0

    if ( callnum == 1  ) then

       if ( MasterProc  ) then
          write(fname,fmt='(''rough.dat'')')
          write(6,*) 'filename for landuse ',fname
       end if
       needed_found=.false.
       call ReadField(IO_ROUGH,fname,r_class,needed_found,fill_needed=.true.)

       ! And convert from real to integer field

       nwp_sea(:,:) = .false.

       if(needed_found)then
          do j=1,ljmax
             do i=1,limax
                if ( nint(r_class(i,j)) == 0 ) nwp_sea(i,j) = .true.
             enddo
          enddo
       endif
    else ! callnum == 2
       if (MasterProc) then
          write(fname,fmt='(''snowc'',i2.2,''.dat'')') current_date%month
          write(6,*) 'filename for snow ',fname

       endif !MasterProc

       needed_found=.not.foundsdepth !needed if not defined through metdata

       call ReadField(IO_SNOW,fname,snow,needed_found)

       if(.not. needed_found)then
             snow=0!should be defined through snowsepth from metdata
       endif
    end if ! callnum == 1

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  end subroutine MetModel_LandUse
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>









  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>







  subroutine tiphys(numt)
    !c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c    written by Trond Iversen,  modified by Hugo Jakobsen, 060994
    !c
    !c    Called from: eulmain.f
    !c
    !c-----------------------------------------------------------------
    !c
    !!        This routine calculates the exner function,
    !!        the  geopotential height, and the vertical exchange coefficient
    !!        in sigma surfaces.
    !!        The height zi of the "well mixed layer" or ABL-height
    !!        is also calculated.
    !c
    !c
    !c    if nroa = 1 also roa is calculated.
    !c    if nfirst=1 for the initial timelevel
    !c
    !c
    !c-----------------------------------------------------------------
    !c    routines called:
    !c
    !c        smoosp
    !c
    !c
    !c-----------------------------------------------------------------
    !c
    !c    DescriPTion of the parameters/variables defined in this file:
    !c
    !c
    !c    absfac    : |xfac|
    !c    abshd    : |fm|
    !c    amax1    : fortran function, choosing largest value
    !c    amin1    : fortran function, choosing smallest value
    !c    CP    : heat capaciyt of air at constant pressure, J/(kg K)
    !c    delq    : available heat flux for developing the unstable ABL, J/m2
    !!              : heat-input per m2 from the ground during unstable BL
    !c    deltaz    : zm(i,k) - zm(i,k+1), m
    !c    dpidth    : heat increasement in accordance with temp. increasement, J/m2
    !c    dth    : iterative increament in potential temperature
    !c    dth0    : accumulated increament in iterative temperature
    !c    dtz    : time interwall for integration of surface heat fluxes
    !c          in the ABL-height calculations, s
    !c    dvdz    : Wind shear, 1/s
    !c    eps    : small number avoiding ri to become infinitely large
    !c    exfrco    : parameter in the Kz model
    !c    exnm    : exner function in the full sigma-levels, J/(kg K)
    !c    exns    : exner function in the half sigma-levels, J/(kg K)
    !c    fh    : surface flux of sensible heat, W/m2
    !c    fl    : surface flux of sensible heat, W/m2 ! ds u7.4vg
    !c    fm    : surface stress (flux of momentum), N/m2
    !c    g    : gravitational acceleration, m/s2
    !c    hs    : height of surface layer (i.e. prandtl-layer), m
    !c    hsl    : (= hs/l, where l is the monin-obhukov length)
    !c    i    : grid index in x-direction
    !c    iip    : limax + 1
    !c    limax    : max number of grid points in x-direction
    !c    iznew    : index for new value of the ABL-height, m
    !c    izold    : index for previous value of the ABL-height, m
    !c    j    : grid index in y-direction
    !c    jjp    : ljmax + 1
    !c    ljmax    : max number of grid points in y-direction
    !c    k    : grid index in vertical-direction
    !c    kkk    : helping index for the cycling of ABL-height
    !c    kkm    : number of full s-levels, *** not used ***
    !c    KMAX_BND    : max number of vertical half levels
    !c          in sigma coordinates
    !c    KMAX_MID    : max number of vertical full levels
    !c          in sigma coordinates
    !c       kzmax   : maximum value of Kz_ms2, m2/s
    !c       kzmin   : minimum value of Kz_ms2, m2/s
    !c    ndth    : do variable for convective ABL-height iteration loop
    !c    nh1    : counts number of layers below zlimax
    !c    nh2    : counts number of layers with Kz > ( Kz )limit
    !c    nr    : number of met.fields stored in arrays (= 1 or 2)
    !c    nt    : time counting variable of the outer time-loop
    !c    p    : local pressure, hPa (mb)
    !c       pi      : pi = 4.*atan(1.) = 3.14 ...
    !c    pidth    : heat used to adjust air temperature, J/m2
    !c    pref    : refference pressure (at ground level), 1.e+5 Pa
    !c    ps    : surface pressure, hPa
    !c    PT    : pressure at the top of the model atmosphere, hPa (mb)
    !c    pz    : local pressure in half sigma levels, hPa (mb),
    !c          helping array (j - slices) for pressure
    !c    pzpbl    : stores H(ABL) for averaging and plotting purposes, m
    !c    ri    : richardson`s number
    !c    ri0    : critical richardson`s number
    !c    risig    : richardson's number in sigmas-levels
    !c    roas    : air density at surface, kg/m3
    !c    sigma_bnd    : height of the half-sigma layers
    !c    sigma_mid    : height of the full-sigma layers
    !c    sm    : height of the surface layer in s-coordinates (4% of H(ABL), m
    !c    th    : potensial temperature (theta), K
    !c    t2_nwp    : potensial temperature at 2m height, K
    !c    thadj    : adjustable surface temperature, K
    !c    thsrf    : potensial temperature at the surface, K
    !c    trc    : helping variable telling whether or not unstable ABL exists
    !!              :       0 => no need for further calc. of ziu
    !!              :       1 => ziu not found yet.
    !c    u_xmj    : wind speed in the x-direction, m/s
    !c    ustar    : friction velocity, m/s
    !c    v_xmi    : wind speed in the y-direction, m/s
    !c       ven     : ventilation coefficient, m3
    !c       venav   : time averaged ventilation coefficient, m3
    !c       venmax  : maximum value of ven, m3
    !c       venmin  : minimum value of ven, m3
    !c       ven00   : averaged ventilation coefficient at 00 UTC, m3
    !c       ven06   : averaged ventilation coefficient at 06 UTC, m3
    !c       ven12   : averaged ventilation coefficient at 12 UTC, m3
    !c       ven18   : averaged ventilation coefficient at 18 UTC, m3
    !c    vdfac    : factor for reduction of vD(1m) to vD(hs)
    !!              : i.e. factor for aerodynamic resistance towards dry deposition
    !!              : vd(50m) = vd(1m)/(1 + vd(1m)*vdfac)
    !c    x12    : mixing length squared, m2
    !c    xfac    : helping variable for reducing concentrations to 1m values
    !c    xfrco    : parameter in the Kz model
    !c    KAPPA    : r/CP (-)
    !c    KARMAN    : von Karmans constant
    !c    xkdz    : the vertical derivative of xkhs at hs, m/s
    !!              : i.e. vertical gradient of xkhs
    !c    xkhs    : diffusivity at hs (in surface layer), m2/s
    !!              : i.e. vertical exchange coeff. on top of prandtl-layer
    !c    Kz_ms2    : estimated exchange coefficient, Kz,  in intermediate
    !c          sigma levels, m2/s
    !c    xksm    : spacially smoothed Kz in z direction, m2/s.
    !!              : Kz_ms2 smoothed over three adjacent layers
    !c    xkzi    : local helping array for the vertical diffusivity, m2/s
    !!              : i.e. vertical exchange coeff. on top of ABL for unstable BL
    !c    zi    : Height of ABL (final value), m
    !c    zlimax    : maximum value of ABL-height, zi, (2000), m
    !c    zimhs    : ziu - hs
    !c    zimin    : minimum value of ABL-height, zi, (200), m
    !c    zimz    : ziu - zs_bnd
    !c    zis    : height of the stable ABL, m
    !c    ziu    : height of the unstable ABL, m
    !c    zixx    : Height og ABL (intermediate value), m
    !c    zm    : geopotential height of full sigma levels above topography, m
    !c    zmhs    : zs_bnd - hs
    !c    zs_bnd    : geopotential height of  half sigma levels above topography, m
    !c    ztop    : height of the uppermost layer in s-coordinates
    !c
    !c-------------------------------------------------------------------
    !c..the following sketches the sigma-surfaces:
    !c
    !c
    !!                ///////////////////
    !c    sigma_bnd(1) = 0 - -sigmas - - - - - sdot(1) = 0, Kz_ms2(1)=xksm(1)=0,
    !!                                           pr(1)=0,PT,exns(1), zs_bnd(1)
    !c
    !c        sigma_mid(1) ---sigmam---------- u_xmj, v_xmi, th, q, cw, exnm (1)
    !c
    !c
    !c        sigma_bnd(2) - - - - s - - - - - sdot(2), Kz_ms2(2), exns(2), pr(2)
    !!                                                 zs_bnd(2), xksm(2)
    !c
    !c        sigma_mid(2) --------m---------- u_xmj, v_xmi, th, q, cw, exnm (2)
    !c
    !c
    !c        sigma_bnd(3) - - - - s - - - - - sdot(3), Kz_ms2(3), exns(3), pr(3)
    !!                                                 zs_bnd(3), xksm(3)
    !c
    !c        sigma_mid(3) --------m---------- u_xmj, v_xmi, th, q, cw, exnm (3)
    !c
    !c
    !c        sigma_bnd(4) - - - - s - - - - - sdot(4), Kz_ms2(4), exns(4), pr(4)
    !!                                                  zs_bnd(4), xksm(4)
    !c
    !c        sigma_mid(4) --------m---------- u_xmj, v_xmi, th, q, cw, exnm (4)
    !c
    !c
    !c        sigma_bnd(5) - - - - s - - - - - sdot(5), Kz_ms2(5), exns(5), pr(5)
    !!                                                  zs_bnd(5), xksm(5)
    !c
    !!                        :
    !!                        :
    !c
    !  sigma_bnd(KMAX_BND-1) - - - - s - - -  - sdot(KMAX_BND-1), Kz_ms2(KMAX_MID),
    !!                                    exns(KMAX_BND-1),zs_bnd(KMAX_BND-1),
    !!                                    pr(KMAX_BND-1),xksm(KMAX_MID)
    !c
    !c    sigma_mid(KMAX_MID) --------m-------- u_xmj, v_xmi, th, q, cw, exnm (KMAX_MID);
    !!                                    this level is assumed to be
    !!                                    the top of Prandtl-layer (LAM50E)
    !c
    ! sigma_bnd(KMAX_BND) = 1- -  - s - - - - sdot(KMAX_BND) = 0, ps, t2_nwp, fh,
    !!                ///////////////////        fm, mslp, Kz_ms2(KMAX_MID)=0,
    !!                                           exns(KMAX_BND), zs_bnd(KMAX_BND),
    !!                                           pr(KMAX_BND),xksm(KMAX_BND)=0.
    !c
    !c
    !c**********************************************************************
    logical, parameter :: DEBUG_KZ = .false.
    logical, parameter :: PIELKE_KZ = .true.    ! Default
    logical, parameter :: TKE_DIFF = .false.  !!! CODE NEEDS TESTING/TIDY UP
    ! STILL!!!!
    !  Kz-tests
    real, parameter :: KZ_MINIMUM = 0.001   ! m2/s
    real, parameter :: KZ_MAXIMUM = 1.0e3 ! m2/s - as old kzmax
    real, parameter :: KZ_SBL_LIMIT = 0.1 ! m2/s - Defines stable BL height

    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID)::exnm
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND)::exns,zs_bnd
    real, dimension(MAXLIMAX,KMAX_MID)::zm,dthdz,deltaz,thc
    real, dimension(MAXLIMAX,KMAX_BND)::risig,xksm,pz
    real, dimension(MAXLIMAX)::zis,delq,thsrf,trc,pidth,dpidth,xkhs,xkdz,xkzi,&
         hs,xkh100
    real,  dimension(MAXLIMAX,MAXLJMAX)::ziu,help,a,zixx,uabs,vdfac
    real ::lim,xdthdz,zmmin,zimin,zlimax,kzmin,kzmax,sm,pref,eps,ric,&
         ric0,dthdzm,dthc,xdth,xfrco,exfrco,hsl,dtz,p,dvdz,xl2,uvhs,zimhs,&
         zimz,zmhs,ux0,fac,fac2,dex12,ro
    real ::h100 ! Top of lowest layer - replaces 100.0
    real ::hsurfl
real ::  p_mid(KMAX_MID), exf2(KMAX_MID), & !TESTzi
         p_bnd(KMAX_BND), & !TESTzi
    ziTI, ziSeibert, ziJericevic, ziVenki, ziZil !TEST

    integer i,j,k,km,km1,km2,kabl,iip,jjp,numt,kp, nr


    integer, dimension(MAXLIMAX) ::  nh1, nh2

    !  Check:
    call CheckStop( KZ_SBL_LIMIT < 1.01*KZ_MINIMUM,   &
         "SBLlimit too low! in Met_ml")

    iip = limax+1
    jjp = ljmax+1

    !
    !     Preliminary definitions
    !
    nr = 2
    if (numt.eq.1) nr = 1
    pref = 1.e+5

    !from ModelC      pi = 4.*atan(1.)

    zlimax = 3000.
    zimin = 100.
    zmmin = 200.

    eps = 0.01
    dtz = 3600.
    sm = 0.04
    !
    !
    !..preset=zero:
    xksm(:,:)   = 0
    risig(:,:)  = 0.
    Kz_m2s(:,:,:)= 0.


    !c..................................
    !c..exner-function in the full sigma-levels..
    !c
    do  k=1,KMAX_MID
       do j=1,ljmax
          do i=1,limax

             !c..pressure (pa)
             p = PT + sigma_mid(k)*(ps(i,j,nr) - PT)
             !c..exner (j/k kg)
             exnm(i,j,k)= CP * Exner_nd(p)
          end do
       end do
    end do

    !c.........................................
    !c..procedure to arrive at mixing height..:
    !c
    !c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !c
    !c     Start j-slice here.
    !c


    do j=1,ljmax

       !c..exner in half-sigma levels:
       do  k=1,KMAX_BND
          do i=1,limax
             p = PT + sigma_bnd(k)*(ps(i,j,nr) - PT)
             pz(i,k) = p
             exns(i,j,k)= CP * Exner_nd(p)
          end do
       end do
       !
       !
       !.. exns(KMAX_BND), th(KMAX_BND) and height of sigmas:
       do  i=1,limax
          zs_bnd(i,j,KMAX_BND)=0.
       end do
       !
       !     Height of the half levels
       !
       do  k=KMAX_BND-1,1,-1
          do i=1,limax
             zs_bnd(i,j,k)=zs_bnd(i,j,k+1)+th(i,j,k,nr)*&
                  (exns(i,j,k+1)-exns(i,j,k))/GRAV

          end do
       end do
       !
       !..height of sigma:
       do  k=1,KMAX_MID
          do i=1,limax
             zm(i,k) = ((exnm(i,j,k)-exns(i,j,k))*zs_bnd(i,j,k+1)&
                  + (exns(i,j,k+1)-exnm(i,j,k))*zs_bnd(i,j,k))&
                  / (exns(i,j,k+1)-exns(i,j,k))
!             zm3d(i,j,k) = zm(i,k)
          end do
       end do

   ! Shows that z_bnd == zs_bnd and z_mid == zm!!!
   !if( DEBUG_HMIX .and. debug_proc .and. & 
   !     !   modulo( current_date%hour, 3)  == 0  .and. &
   !        current_date%seconds == 0 .and. j == debug_jloc ) then
   !   i = debug_iloc !,    j = debug_jloc
   !    do k = 17, 1, -4
   !     write(6,"(a,i3,4f15.4)") "DEBUG_Hmix Zs", k, z_bnd(i,j,k), zs_bnd(i,j,k), z_mid(i,j,k), zm(i,k)
   !    end do
   ! end if


       !----------------------------------------------------------------------
       !...........................................
       !..the following variables in sigmas-levels:
       !
       do  k=2,KMAX_MID
          km=k-1
          do i=1,limax
             !
             !.........................
             !..wind sheare
             !
             ! Slightly different formulation of dvdz than in metvar

             !BUG dvdz = ( (u_xmj(i,j,km,nr)-u_xmj(i,j,k,nr))**2 &
             !BUG      + (v_xmi(i,j,km,nr)-v_xmi(i,j,k,nr))**2 + eps)
              dvdz = (u_mid(i,j,km)-u_mid(i,j,k))**2 + & ! Actually (dv)^2 
                     (v_mid(i,j,km)-v_mid(i,j,k))**2 + eps  
             !
             risig(i,k)=(2.*GRAV/(th(i,j,km,nr)+th(i,j,k,nr)))*&
                  (th(i,j,km,nr)-th(i,j,k,nr))*(zm(i,km)-zm(i,k))&
                  /dvdz
             !........................
             !..mixing length squared:
             !
             xl2=(KARMAN*amin1(zs_bnd(i,j,k),zmmin))**2

             !
             !..............................
             !..critical richardsons number:
             !
             ric0=0.115*((zm(i,km)-zm(i,k))*100.)**0.175
             ric=amax1(0.25,ric0)



             dvdz = sqrt(dvdz)/(zm(i,km)-zm(i,k))

             !..................................................................
             !..exchange coefficient (Pielke,...)
             if ( PIELKE_KZ  ) then
                if (risig(i,k) > ric ) then
                   Kz_m2s(i,j,k) = KZ_MINIMUM
                else
                   Kz_m2s(i,j,k) = 1.1 * (ric-risig(i,k)) * xl2 * dvdz /ric
                end if
             else

                !..exchange coefficient (blackadar, 1979; iversen & nordeng, 1987):
                !
                if(risig(i,k).le.0.) then
                   Kz_m2s(i,j,k)=xl2*dvdz*sqrt(1.1-87.*risig(i,k))
                elseif(risig(i,k).le.0.5*ric) then
                   Kz_m2s(i,j,k)=xl2*dvdz*(1.1-1.2*risig(i,k)/ric)
                elseif(risig(i,k).le.ric) then
                   Kz_m2s(i,j,k)=xl2*dvdz*(1.-risig(i,k)/ric)
                else
                   Kz_m2s(i,j,k)=0.001
                endif
             end if ! Pielke or Blackadar
             !
          end do
       end do
       !
       !tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
       !
       !---------------------------------------------------------------------
       !..................................
       !..height of stable boundary layer:
       !
       !.........................................................
       !..vertical smoothing of Kz_ms2 over three adjacent layers:
       !
       k=2
       km=1
       kp=3
       do i=1,limax
          xksm(i,k)=( (zm(i,km)-zm(i,k))*Kz_m2s(i,j,k)&
               + (zm(i,k)-zm(i,kp))*Kz_m2s(i,j,kp) )&
               / ( zm(i,km) - zm(i,kp) )
       enddo
       !c
       k=KMAX_MID
       km2=k-2
       km1=k-1
       do i=1,limax
          xksm(i,k)=( (zm(i,km2)-zm(i,km1))*Kz_m2s(i,j,km1)&
               + (zm(i,km1)-zm(i,k))*Kz_m2s(i,j,k) )&
               / ( zm(i,km2) - zm(i,k) )
       enddo
       !c
       do k = 3,KMAX_MID-1
          km1=k-1
          km2=k-2
          kp=k+1
          do i=1,limax
             xksm(i,k)=(  (zm(i,km2)-zm(i,km1))*Kz_m2s(i,j,km1)&
                  + (zm(i,km1)-zm(i,k))*Kz_m2s(i,j,k)&
                  + (zm(i,k)-zm(i,kp))*Kz_m2s(i,j,kp) )&
                  / ( zm(i,km2) - zm(i,kp) )
          enddo
       enddo


       !c
       !c............................................................
       !c..The height of the stable BL is the lowest level for which:
       !c..xksm .le. 1 m2/s (this limit may be changed):
       !c
       do i = 1,limax
          zis(i)=zimin
          nh1(i) = KMAX_MID
          nh2(i) = 1
       enddo
       !c
       do k=KMAX_MID,2,-1
          do i=1,limax

             if(xksm(i,k) >= KZ_SBL_LIMIT .and. nh2(i) == 1) then
                nh1(i)=k   ! Still unstable
             else
                nh2(i)=0   ! Now stable
             endif
             !c
          end do
       end do
       !c
       do i=1,limax
          !c
          k=nh1(i)
          !c
          if(zs_bnd(i,j,nh1(i)).ge.zimin) then

             if( abs(xksm(i,k)-xksm(i,k-1)) .gt. eps) then

                zis(i)=((xksm(i,k)-KZ_SBL_LIMIT )*zs_bnd(i,j,k-1) &
                     + (KZ_SBL_LIMIT -xksm(i,k-1))*zs_bnd(i,j,k))&
                     /(xksm(i,k)-xksm(i,k-1))
             else

                zis(i)=zimin
             endif

          endif
          !
       end do
       !
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !
       !---------------------------------------------------------------------
       !....................................
       !..height of unstable boundary layer:
       !
       !
       !..assuring that th is increasing with height.
       !..adjusted th-sounding is assigned to thc-array.
       !..This adjusted th is not meant to be used in
       !..other parts of the model program
       !
       dthdzm = 1.e-4
       do i =1,limax
          thc(i,KMAX_MID)=th(i,j,KMAX_MID,nr)
          do k=KMAX_MID-1,1,-1

             dthc = (th(i,j,k,nr)-th(i,j,k+1,nr))&
                  / (zm(i,k)-zm(i,k+1))

             dthdz(i,k)=amax1(dthc,dthdzm)

             thc(i,k)=thc(i,k+1)+dthdz(i,k)*(zm(i,k)-zm(i,k+1))

          enddo
       enddo

       !
       !
       !..estimated as the height to which an hour's input
       !..of heat from the ground is vertically distributed,
       !..assuming dry adiabatic adjustment.
       !
       !
       do  i=1,limax

          delq(i)=-amin1((fh(i,j,nr)),0.)*dtz
          thsrf(i)=0.
          ziu(i,j)=0.
          !
          !.................................
          !..trc=1 for unstable BL (delq>0):
          !..   =0 for stable BL (delq=0):
          !
          if(delq(i).gt.0.00001) then
             trc(i)=1.
          else
             trc(i)=0.
          endif
          !
          !------------------------------------------------------------
          ! calculating the height of unstable ABL
          !
          !!      if(trc(i).eq.1.) then

          kabl = KMAX_MID

          do while( trc(i).eq.1)
             !! 28     if(trc(i).eq.1.) then
             kabl = kabl-1
             pidth(i)=0.

             do k=KMAX_MID,kabl,-1
                xdth = thc(i,kabl)-thc(i,k)
                dpidth(i) = exnm(i,j,k)*xdth*(pz(i,k+1)-pz(i,k))/GRAV
                pidth(i) = pidth(i) + dpidth(i)

             end do


         if(pidth(i).ge.delq(i).and.trc(i).eq.1.) then

                !c    at level kabl or below level kabl and above level kabl+1


                thsrf(i) = thc(i,kabl)     &
                     - (thc(i,kabl)-thc(i,KMAX_MID))    &
                     * (pidth(i)-delq(i))/pidth(i)

                xdthdz = (thc(i,kabl)-thc(i,kabl+1))        &
                     / (zm(i,kabl)-zm(i,kabl+1))

                ziu(i,j) = zm(i,kabl+1)                     &
                     + (thsrf(i)-thc(i,kabl+1))/xdthdz

                trc(i)=0.

             endif


             if(kabl.le.4 .and. trc(i).eq.1.) then

                write(6,*)'Metml ziu calculations failed!'

                ziu(i,j)=zlimax

                trc(i)=0.
             endif

          end do ! while

       end do

       !..iteration, finding height of unstable BL  finished
       !.....................................................................


       do i=1,limax

      if( i_fdom(i) == DEBUG_I .and. j_fdom(j) == DEBUG_J &
            .and.  debug_proc == .true. ) &
           write(6,"(a,3f9.3)") "TIDEBUG zius ", ziu(i,j),zis(i), zixx(i,j)
          zixx(i,j)=amax1(ziu(i,j),zis(i))
          zixx(i,j)=amin1(zlimax,zixx(i,j))
       end do



    end do !     End j-slice
    !!-------------------------------------------




    !..spatial smoothing of new zi:

    call smoosp(zixx,zimin,zlimax)

    do j=1,ljmax
       do i=1,limax
          pzpbl(i,j) = zixx(i,j)
       enddo
    enddo

   if( DEBUG_HMIX .and. debug_proc .and. & 
         modulo( current_date%hour, 3)  == 0  .and. &
           current_date%seconds == 0  ) then

      i = debug_iloc
      j = debug_jloc

      call SeibertRiB_Hmix(&
          u_mid(i,j,KPBL_MAX:KMAX_MID),  &
          v_mid(i,j,KPBL_MAX:KMAX_MID),  &
          z_mid(i,j,:),  &
          th(i,j,:,nr),  &
          ziSeibert)
      call JericevicRiB_Hmix(&
          u_mid(i,j,KPBL_MAX:KMAX_MID),  &
          v_mid(i,j,KPBL_MAX:KMAX_MID),  &
          z_mid(i,j,:),  &
          th(i,j,:,nr),  &
          ziJericevic)
      ziVenki = Venkatram_Hmix(ustar_nwp(i,j))
      ! no 1/L yet
      !call Zilitinkevich_Hmix(&
      !    ustar_nwp(i,j),  &
      !    ustar_nwp(i,j),  &
      !    ziVenki)

     !   "old" exner-function of the full-levels

      p_bnd(:) = sigma_bnd(:)*(ps(i,j,nr) - PT) + PT
      p_mid(:) = sigma_mid(:)*(ps(i,j,nr) - PT) + PT
      exf2(:) = CP * Exner_nd(p_mid(:))
      call TI_BLphysics(&
          u_mid(i,j,KPBL_MAX:KMAX_MID),  &
          v_mid(i,j,KPBL_MAX:KMAX_MID),  &
          z_mid(i,j,:),  &
          z_bnd(i,j,:),  &
          fh(i,j,nr),  &
          th(i,j,:,nr),  &
          exf2(:),  &
          p_bnd(:), &
          ziTI,.true.)

      write(6,"(a,4i5, f8.2,11f7.1)") "DEBUG_Hmix: ", &
           current_date%month,&
           current_date%day,&
           current_date%hour,&
           current_date%seconds, fh(i,j,nr), &
           ziu(i,j),zis(i),  &
           pzpbl(i,j), ziTI, ziSeibert, ziJericevic, ziVenki
    end if


    !cttttttttttttttttttttttttttttttttttttttttttttttttttttttt
    !c..height of ABL finished..............................
    !c------------------------------------------------------



    !----------------------------------------------------------!
    if( TKE_DIFF ) then
       call tkediff (nr)                            ! guta
    else if (NWP_Kz) then
       do  k=2,KMAX_MID
         do i=1,limax
            do j=1,ljmax
! hb KZ_met are accumulated values over three hours, to get m2/s divide with 60*60*3  
               SigmaKz(i,j,k,nr)=Kz_met(i,j,k,nr)/(60*60*3)
               Kz_m2s(i,j,k)=SigmaKz(i,j,k,nr)*SigmaKz_2_m2s(roa(i,j,k,nr),ps(i,j,nr))
            enddo
         enddo
       enddo

     else if (OBRIAN_Kz ) then

       call O_Brian(nr, KZ_MINIMUM, KZ_MAXIMUM, zimin, zs_bnd, ziu  &
            , exns, exnm, zixx, DEBUG_HMIX )
    end if

   if( MasterProc .and. DEBUG_Kz)then            
       write(6,*)               &
       '*** After Set SigmaKz', sum(SigmaKz(:,:,:,nr)), minval(SigmaKz(:,:,:,nr)), maxval(SigmaKz(:,:,:,nr)), DEBUG_Kz, 'NWP_Kz:',NWP_Kz,'OBRIAN_Kz:',OBRIAN_Kz, '*** After convert to z',sum(Kz_m2s(:,:,:)), minval(Kz_m2s(:,:,:)), maxval(Kz_m2s(:,:,:))
    endif

    !----------------------------------------------------------!
  if( DEBUG_Kz .and. debug_proc .and. &
         modulo( current_date%hour, 3)  == 0  .and. &
           current_date%seconds == 0  ) then
    
      i = debug_iloc
      j = debug_jloc
      do k=1,KMAX_MID

      write(6,"(a,4i5, 4f17.12)") "DEBUG_Kz: ", &
           current_date%month,&
           current_date%day,&
           current_date%hour,&
           k,&
           SigmaKz(i,j,k,nr),&
           Kz_m2s(i,j,k)
      enddo
   end if

  end subroutine tiphys

  !c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !





  subroutine smoosp(f,rmin,rmax)

    !c    file: eulmet-mnd.f
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c    written by Trond Iversen,  modified by Hugo Jakobsen, 080994
    !       parallellized and modified by Peter February 2003
    !
    !c
    !c    Called from: tiphys.f
    !c
    !c----------------------------------------------------------------------
    !c
    !c    This routine applies the shapiro filter with s=0.5 and s=-0.5
    !c    to the field f usinh h as a work space also the boundaries
    !c    are smoothed. f contains the smoothed field upon return.
    !c

    !c    Definition of the variables:
    !c
    !c
    !c    f    : data to be smoothed
    !c    iif    : =limax
    !c    jjf    : =ljmax
    !c    h1,h2    : = help variable
    !c    rmin    : min allowed
    !c    rmax    : max allowed
    !c
    implicit none

    real, intent(inout) :: f(MAXLIMAX,MAXLJMAX)
    real, intent(in)    :: rmin,rmax

    real, dimension(MAXLIMAX+4,MAXLJMAX+4) :: h1, h2
    real, dimension(MAXLIMAX,2)            :: f_south,f_north
    real, dimension(MAXLJMAX+2*2,2)        :: f_west,f_east
    real s

    integer  thick
    integer iif,jjf,is,i,j,ii,jj,iifl,jjfl

    iif=limax
    jjf=ljmax

    thick=2  !we fetch 2 neighbors at once, so that we don't need to call
    ! readneighbours twice
    iifl=iif+2*thick
    jjfl=jjf+2*thick

    call readneighbors(f,f_south,f_north,f_west,f_east,thick)

    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f(i,j)
       enddo
    enddo
    do j=1,thick
       do i=1,iif
          ii=i+thick
          h1(ii,j) = f_south(i,j)
       enddo
    enddo

    do j=1,thick
       jj=j+jjf+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f_north(i,j)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          h1(i,j) = f_west(j,i)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          ii=i+iif+thick
          h1(ii,j) = f_east(j,i)
       enddo
    enddo

    do j=1,jjfl
       h2(1,j) = 0.
       h2(iifl,j) = 0.
    enddo

    do i=1,iifl
       h2(i,1) = 0.
       h2(i,jjfl) = 0.
    enddo
    !! 44 format(I2,30F5.0)

    do  is=2,1,-1

       s=is-1.5  !s=0,5 s=-0.5
       if(is /= 2)h1=h2

       !..the smoothing

       do  j=2,jjfl-1
          do  i=2,iifl-1
             h2(i,j)=(1.-2.*s+s*s)*h1(i,j)&
                  + 0.5*s*(1.-s)*(h1(i+1,j)+h1(i-1,j)+h1(i,j+1)+h1(i,j-1))  &
                  + s*s*(h1(i+1,j+1)+h1(i-1,j-1)+h1(i+1,j-1)+h1(i-1,j+1))/4.
             h2(i,j) = amax1(h2(i,j),rmin)
             h2(i,j) = amin1(h2(i,j),rmax)
          end do
       end do

    end do


    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          f(i,j)=h2(ii,jj)
       enddo
    enddo

  end subroutine smoosp
  !  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




  subroutine readneighbors(data,data_south,data_north,data_west,data_east,thick)


    ! Read data at the other side of the boundaries
    !
    ! thick is the number of gridcells in each direction to be transferred
    ! Note that we also fetch data from processors in the "diagonal"
    ! directions
    !
    ! Written by Peter February 2003
    !
    !Note,
    !The data_west(jj,:)=data(1,j) is not a bug: when there is no west
    !neighbour,
    !the data is simply copied from the nearest points: data_west(jj,:) should
    !be =data(-thick+1:0,j), but since this data does not exist, we
    !put it =data(1,j).


    implicit none

    integer, intent(in) :: thick
    real,intent(in), dimension(MAXLIMAX,MAXLJMAX) ::data
    real,intent(out), dimension(MAXLIMAX,thick) ::data_south,data_north
    real,intent(out), dimension(MAXLJMAX+2*thick,thick) ::data_west,data_east

    integer :: msgnr,info
    integer :: i,j,tj,jj,jt

    !check that limax and ljmax are large enough
    call CheckStop(limax < thick, "ERROR readneighbors in Met_ml")
    call CheckStop(ljmax < thick, "ERROR readneighbors in Met_ml")


    msgnr=1

    data_south(:,:)=data(:,1:thick)
    data_north(:,:)=data(:,ljmax-thick+1:ljmax)
    if(neighbor(SOUTH) >= 0 )then
       CALL MPI_SEND( data_south , 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(SOUTH), msgnr, MPI_COMM_WORLD, INFO)
    endif
    if(neighbor(NORTH) >= 0 )then
       CALL MPI_SEND( data_north , 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(NORTH), msgnr+9, MPI_COMM_WORLD, INFO)
    endif

    if(neighbor(SOUTH) >= 0 )then
       CALL MPI_RECV( data_south, 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(SOUTH), msgnr+9, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       do tj=1,thick
          data_south(:,tj)=data(:,1)
       enddo
    endif
44  format(I2,30F5.0)
    if(neighbor(NORTH) >= 0 )then
       CALL MPI_RECV( data_north, 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(NORTH), msgnr, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       do tj=1,thick
          data_north(:,tj)=data(:,ljmax)
       enddo
    endif

    jj=0
    do jt=1,thick
       jj=jj+1
       data_west(jj,:)=data_south(1:thick,jt)
       data_east(jj,:)=data_south(limax-thick+1:limax,jt)
    enddo
    do j=1,ljmax
       jj=jj+1
       data_west(jj,:)=data(1:thick,j)
       data_east(jj,:)=data(limax-thick+1:limax,j)
    enddo
    do jt=1,thick
       jj=jj+1
       data_west(jj,:)=data_north(1:thick,jt)
       data_east(jj,:)=data_north(limax-thick+1:limax,jt)
    enddo

    if(neighbor(WEST) >= 0 )then
       CALL MPI_SEND( data_west , 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(WEST), msgnr+3, MPI_COMM_WORLD, INFO)
    endif
    if(neighbor(EAST) >= 0 )then
       CALL MPI_SEND( data_east , 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(EAST), msgnr+7, MPI_COMM_WORLD, INFO)
    endif

    if(neighbor(WEST) >= 0 )then
       CALL MPI_RECV( data_west, 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(WEST), msgnr+7, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       jj=0
       do jt=1,thick
          jj=jj+1
          data_west(jj,:)=data_south(1,jt)
       enddo
       do j=1,ljmax
          jj=jj+1
          data_west(jj,:)=data(1,j)
       enddo
       do jt=1,thick
          jj=jj+1
          data_west(jj,:)=data_north(1,jt)
       enddo
    endif
    if(neighbor(EAST) >= 0 )then
       CALL MPI_RECV( data_east, 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE, &
            neighbor(EAST), msgnr+3, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       jj=0
       do jt=1,thick
          jj=jj+1
          data_east(jj,:)=data_south(limax,jt)
       enddo
       do j=1,ljmax
          jj=jj+1
          data_east(jj,:)=data(limax,j)
       enddo
       do jt=1,thick
          jj=jj+1
          data_east(jj,:)=data_north(limax,jt)
       enddo
    endif

  end subroutine readneighbors






  !************************************************************************!
  subroutine tkediff (nr)                                             !
    !************************************************************************!
    !                                                                        !
    !    This routine computes vertical eddy diffusivities as a function     !
    !    altitude, height of PBL, and a velocity scale, square root of       !
    !    turbulent kinetic energy (TKE). This is a non-local scheme.         !
    !    The TKE at the surface is diagnosed using scales for horizontaland  !
    !    vertical velocities (ustar and wstar) in the surface layer          !
    !    (Alapaty 2004; Holstag et al. 1990 and Mihailovic et al. 2004)      !
    !    PBL ht is calculated using the EMEP formulation                     !
    !                                                                        !
    !    Written by DT Mihailovic (October 2004)                             !
    !    EMEP polishing and comments: JE Jonson and P Wind                   !
    !************************************************************************!

    implicit none

    !     Local constants
    real   , parameter :: SZKM=1600.     &   ! Constant (Blackadar, 1976)
         ,CKZ=0.001      &   ! Constant (Zhang and Athens, 1982)
         ,REFPR=1.0E+05  &   ! Referent pressure
         ,KZ0LT=1.0E-04  &   ! Constant (Alapaty et al., 1997)
         ,RIC=0.10       &   ! Critical Richardson number
                                ! (Holstlag et al., 1993)
         ,ROVG=RGAS_KG/GRAV        ! Used in Calculation of R-number
    integer, parameter :: KLM =KMAX_MID-1



    !     INPUT
    integer, intent(in) :: nr  ! Number of meteorological stored
    ! in arrays (1 or 2)

    !     OUTPUT
    !     skh(i,j,k,nr) array
    !     Values of the Kz coefficients (eddyz (i,j,k)) are transformed nto
    !     sigma system and then they stored in this array which is later used
    !     in ADVECTION module


    !     Local arrays

    integer, dimension(MAXLIMAX,MAXLJMAX)      :: iblht   ! Level of the PBL top
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND):: eddyz   ! Eddy coefficients
    ! (m2/s)
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID):: &
         t_virt    &! Potential temperature (K)
         ,e         &! Kinetic energy with respect to height (m2/s2)
         ,dzq       &! Thickness of sigma interface layers (m)
         ,u_mid     &! Wind speed in x-direction (m/s)
         ,v_mid      ! Wind speed in y-direction (m/s)

    real, dimension(MAXLIMAX,MAXLJMAX,KLM):: &
         dza         ! Thickness of half sigma layers (m)

    real, dimension(MAXLIMAX,MAXLJMAX):: &
         pblht ,    &! PBL (Holstag, 1990) (m)
         h_flux,    &! Sensible heat flux  (W/m2)
         ust_r ,    &! Friction velocity (m/s)
         mol   ,    &! Monin-obukhov length (m)
         wstar       ! Convective velocity (m/s)

    real, dimension(KMAX_BND) :: rib         ! Bulk Richardson number

    real, dimension(KMAX_MID) :: &
         rich,       &! Richardson number
         psi_zi       ! Used in the vertical integration

    real, dimension (10) ::      &
         psi_z       & ! Used for calculating
         , zovh          ! TKE

    !     Local variables
    real dtmp, tog, wssq1, wssq2, wssq, tconv, wss, wst, PSI_TKE,    &
         dusq, dvsq, ri, ss, dthdz, busfc, zvh,                   &
         part1, part2, fract1, fract2, apbl, fac0, fac02, kz0,    &
         cell, dum1, rpsb, press, teta_h, u_s, goth, pressure

    integer i, j, k, l, kcbl


    !     Functions for averaging the vertical turbulent kinetic energy
    !      (Alapaty, 2003)
    data psi_z /0.00,2.00,1.85,1.51,1.48,1.52,1.43,1.10,1.20,0.25/
    data zovh  /0.00,0.05,0.10,0.20,0.40,0.60,0.80,1.00,1.10,1.20/

    !     Store the NMW meteorology and variables derived from its

    !     Change the sign
    h_flux(1:limax,1:ljmax)=-fh(1:limax,1:ljmax,nr)

    !     Avoid devision by zero later in the code

    where (ABS(h_flux(1:limax,1:ljmax))<0.0001) h_flux(1:limax,1:ljmax)=0.0001

    !     Check PBL height   ! strange tests! Negative pzpbl check? From 1 to 100m
    !   - odd!
    do i=1,limax
       do j=1,ljmax
          if(ABS(pzpbl(i,j)) < 1.) then
             pzpbl(i,j)=100.
          endif
       enddo
    enddo

    !     Calculate velocity components in the (h) poits (Arakawa notation)
    do k=1,KMAX_MID
       do i=1,limax
          do j=1,ljmax
             !          u_mid(i,j,k)=0.5*(u_xmj(i-1,j  ,k,nr)+u_xmj(i,j,k,nr))
             !          v_mid(i,j,k)=0.5*(v_xmi(i  ,j-1,k,nr)+v_xmi(i,j,k,nr))

             u_mid(i,j,k)=u_xmj(i,j  ,k,nr)
             v_mid(i,j,k)=v_xmi(i  ,j,k,nr)

          enddo
       enddo
    enddo

    !     Avoid small values
    where (ABS(u_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
         u_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001
    where (ABS(v_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
         v_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001

    !     Initialize eddy difussivity arrays
    eddyz(1:limax,1:ljmax,1:KMAX_MID)=0.

    !     Calculate tickness of the full layers
    dzq(1:limax,1:ljmax,1:KMAX_MID) = z_bnd(1:limax,1:ljmax,1:KMAX_MID)  &
         - z_bnd(1:limax,1:ljmax,2:KMAX_BND)

    !     ... and the half sigma layers
    dza(1:limax,1:ljmax,1:KLM) = z_mid(1:limax,1:ljmax,1:KLM)          &
         - z_mid(1:limax,1:ljmax,2:KMAX_MID)

    !     Calculate virtual temperature

    t_virt(1:limax,1:ljmax,1:KMAX_MID) = th(1:limax,1:ljmax,1:KMAX_MID,nr)  &
         * (1.0+0.622*q(1:limax,1:ljmax,1:KMAX_MID,nr))


    !     Calculate Monin-Obuhkov length   (Garratt, 1994)

    do i=1,limax
       do j=1,ljmax
          u_s = ustar_nwp(i,j)
          mol(i,j) = -(ps(i,j,nr)*u_s*u_s*u_s)/                        &
               (KARMAN*GRAV*h_flux(i,j)*KAPPA)
       enddo
    enddo

    !     Calculate the convective velocity (wstar)
    do i=1,limax
       do j=1,ljmax
          wstar(i,j) = GRAV*h_flux(i,j)*pzpbl(i,j)/rho_surf(i,j)    &
               /CP/th(i,j,KMAX_MID,nr)
          if(wstar(i,j) < 0.) then
             wstar(i,j)=-ABS(wstar(i,j))**(0.3333)
          else
             wstar(i,j)=(wstar(i,j))**(0.3333)
          endif
       enddo
    enddo

    !                            ------------------------------------------>
    !     Start with a long loop ------------------------------------------>
    !                            ------------------------------------------>
    DO  i=1,limax
       DO  j=1,ljmax

          rib(1:KMAX_MID) = 0.0        ! Initialize bulk Richardson number

          part1=ust_r(i,j)*ust_r(i,j)*ust_r(i,j)
          wst=AMAX1(wstar(i,j),1.0E-20)
          part2=0.6*wst*wst*wst
          wss=AMAX1(1.0E-4,(part1+part2))
          wss=EXP(0.333333*ALOG(wss))

          if (h_flux(i,j) < 0.0) then
             tconv=0.0                                   ! Holstlag et al. (1990)
          else
             tconv=8.5*h_flux(i,j)/rho_surf(i,j)/CP/wss   ! Conversion to
             ! kinematic flux
          endif

          do k=KMAX_MID,1,-1
             dtmp=t_virt(i,j,k)-t_virt(i,j,KMAX_MID)-tconv
             tog=0.5*(t_virt(i,j,k)+t_virt(i,j,KMAX_MID))/GRAV
             wssq1=u_mid(i,j,k)*u_mid(i,j,k)
             wssq2=v_mid(i,j,k)*v_mid(i,j,k)
             wssq=wssq1+wssq2
             wssq=AMAX1(wssq,1.0E-4)
             rib(k)=z_mid(i,j,k)*dtmp/(tog*wssq)
             if(rib(k).ge.RIC) go to 9001
          enddo
9001      continue

          !     Calculate PBL height according to Holtslag et al. (1993)
          pblht(i,j)=0.
          if(k.ne.KMAX_MID) then
             fract1=(RIC-rib(k+1))/(rib(k)-rib(k+1))
             fract2=1.-fract1
             apbl=z_mid(i,j,k)*fract1
             pblht(i,j)=apbl+z_mid(i,j,k+1)*fract2
             if(pblht(i,j) > z_bnd(i,j,k+1)) then
                kcbl=k
             else
                kcbl=k+1
             endif
          endif
          iblht(i,j)=kcbl

          if(pblht(i,j)<z_bnd(i,j,KMAX_MID)) then
             pblht(i,j)=z_bnd(i,j,KMAX_MID)
             iblht(i,j)=KMAX_MID
          endif


          if(pblht(i,j).le.100.) then              !Minimum of PBL height
             pblht(i,j)=100.
          endif

          !     Find the critical Richardson number (Shir and Borestein, 1976)
          do k=2,iblht(i,j)-1
             rich(k)=0.257*dza(i,j,k)**0.175
          enddo

          !     Free troposphere and cloudy case Kz values estimation
          do k = 2,iblht(i,j)-1
             dusq = (u_mid(i,j,k-1)-u_mid(i,j,k))*(u_mid(i,j,k-1)-u_mid(i,j,k))
             dvsq = (v_mid(i,j,k-1)-v_mid(i,j,k))*(v_mid(i,j,k-1)-v_mid(i,j,k))
             ss = (dusq+dvsq)/(dza(i,j,k-1)*dza(i,j,k-1))+1.E-9
             goth = 2.*GRAV/(t_virt(i,j,k-1)+t_virt(i,j,k))
             dthdz = (t_virt(i,j,k-1)-t_virt(i,j,k))/dza(i,j,k-1)
             ri = goth*dthdz/ss

             !     (Duran and Clemp, 1982)

             kz0 = CKZ*dzq(i,j,k)
             if (ri-rich(k) > 0.) then
                eddyz(i,j,k)=kz0
             else
                eddyz(i,j,k)=kz0+SZKM*SQRT(ss)*(rich(k)-ri)/rich(k)
             endif
             eddyz(i,j,k)=AMIN1(eddyz(i,j,k),100.)
          enddo

          !     Eddy diffusivity coefficients for all regimes in the mixed layer

          do  k=iblht(i,j),KMAX_MID
             if (mol(i,j) < 0.0) then                 !Unstable conditions
                ri=(1.0-15.*z_mid(i,j,k)/mol(i,j))**(-0.25)
                ri=ri/KARMAN/z_mid(i,j,k)
                ri=ri*AMAX1(0.0,pblht(i,j)-z_mid(i,j,k))
                dthdz=ri*ust_r(i,j)**3.
                goth=AMAX1(wstar(i,j),0.0)
                dusq=0.4*goth**3.
                ri=(dthdz+dusq)**(2./3.)
                e(i,j,k)=0.5*ri*(2.6)**(2./3.)        !Moeng and Sullivan (1994)
             else
                ri=z_bnd(i,j,k)/pblht(i,j)               !Stable
                ri=z_mid(i,j,k)/pblht(i,j)               !New
                ri=(1.0-ri)
                ri=AMAX1(0.0,ri)
                ri=(1.0-ri)**1.75
                e(i,j,k)=6.*ust_r(i,j)*ust_r(i,j)*ri  !Lenshow(1988)
             endif

             !     Calculate Ksi function using interpolation in the vertical
             !     Alapaty (2001, 2003)

             zvh=z_mid(i,j,k)/pblht(i,j)
             do l=1,9
                if (zvh > zovh(l).and. zvh < zovh(l+1)) then
                   psi_zi(k)=(psi_z(l+1)-psi_z(l))/(zovh(l+1)-zovh(l))
                   psi_zi(k)=psi_zi(k)*(zvh-zovh(l))
                   psi_zi(k)=psi_zi(k)+psi_z(l)
                   psi_zi(k)=psi_zi(k)/2.0               !Normalized the value
                endif
             enddo
          enddo

          !      Calculate integral for Ksi
          psi_tke=0.
          do k=KMAX_MID,iblht(i,j),-1
             psi_tke=psi_tke+psi_zi(k)*dzq(i,j,k)*sqrt(e(i,j,k))
          enddo

          psi_tke=psi_tke/pblht(i,j)



          do k=iblht(i,j),KMAX_MID          !Calculate coefficients
             goth=psi_tke
             goth=goth*KARMAN*z_mid(i,j,k)
             dthdz=z_mid(i,j,k)/pblht(i,j)
             dthdz=1.0-dthdz
             dthdz=AMAX1(1.0E-2,dthdz)
             if(mol(i,j) > 0.0) then                         !Stable
                goth=sqrt(e(i,j,iblht(i,j)))                 ! Mihailovic (2004)
                goth=goth*KARMAN*z_mid(i,j,k)                 ! -----------------
                dthdz=z_mid(i,j,k)/pzpbl(i,j)                 ! -----------------
                dthdz=1.0-dthdz
                dthdz=AMAX1(1.0E-2,dthdz)
                busfc=0.74+4.7*z_mid(i,j,KMAX_MID)/mol(i,j)
                busfc=AMAX1(busfc,1.0)
                dthdz=dthdz**1.50                                  !test (2004)
                eddyz(i,j,k)=goth*dthdz/busfc
             else
                dthdz=dthdz*dthdz
                busfc=1.0
                eddyz(i,j,k)=goth*dthdz/busfc
             endif
          enddo

          !      Checking procedure
          do k=2,iblht(i,j)-1
             if(eddyz(i,j,k).le.0.0) THEN
                eddyz(i,j,k)= KZ0LT
             endif
          enddo

          !      Avoid phisically unrealistic values
          do k=2,KMAX_MID
             IF(eddyz(i,j,k).le.0.1) then
                eddyz(i,j,k)=0.1
             endif
          enddo

          !     To avoid loss of mass/energy through top of the model
          !     put eddyz (I,J,K) to zero at the last  level from top
          eddyz(i,j,KMAX_BND)=0.0

          !     Calculate eddy coefficients at the interfaces
          do k=2,KMAX_MID
             eddyz(i,j,k)=0.5*(eddyz(i,j,k-1)+eddyz(i,j,k)) !!

             !              if(i.eq.10.and.j.eq.10.) then
             !             if (abs(u_xmj(i,j  ,k,nr)-u_mid(i,j,k)).gt.5.) then
             !
             !       print *,"NEW ",i,j,u_xmj(i,j  ,KMAX_MID,nr),u_mid(i,j,KMAX_MID)
             !              endif
          enddo

          !     Transform values of the eddy coeficients into the the sigma coordinate

          do k=2,KMAX_MID
             eddyz(i,j,k)=eddyz(i,j,k)*((sigma_mid(k)-sigma_mid(    k-1))/   &
                  (    z_mid(i,j,k)-z_mid(i,j,k-1)))**2.

          enddo

       ENDDO                  !---------------------------------------->
    ENDDO                  !---------------------------------------->
    !---------------------------------------->

    !     Store diffusivity coefficients into skh(i,j,k,nr) array
    do k=2,KMAX_MID
       do i=1,limax
          do j=1,ljmax
             SigmaKz(i,j,k,nr)=eddyz(i,j,k)
          enddo
       enddo
    enddo

    ! For plotting set pblht  =  pzpbl

    pzpbl(:,:) = pblht(:,:)

  end subroutine tkediff
  !---------------------------------------------------------------





  !************************************************************************!
  subroutine O_Brian(nr, KZ_MINIMUM, KZ_MAXIMUM, zimin, zs_bnd, ziu  &
       , exns, exnm, zixx, debug_flag )                            !
    !************************************************************************!

    !......................................................
    !..exchange coefficients for convective boundary layer:
    !..o'brien's profile formula:
    !..and the air density at ground level:
    !
    !..constants for free-convection limit:
    !

    integer, intent(in) :: nr

    real, intent(in) :: zimin        &
         ,KZ_MINIMUM   &
         ,KZ_MAXIMUM

    real,intent(in), dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: zs_bnd &
         ,exns
    real,intent(in), dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: exnm

    real,intent(in), dimension(MAXLIMAX,MAXLJMAX) :: ziu    &
         ,zixx
    logical, intent(in) :: debug_flag

    real ::  KBW, invL

    real :: h100   & ! Top of lowest layer - replaces 100.0
         ,xfrco  &
         ,exfrco &
         ,sm     &
         ,ux0    &   ! local ustar
         ,ux3    &   ! ustar**3, ds apr2005
         ,hsl    &
         ,hsurfl &
         ,zimhs  &
         ,zimz   &
         ,zmhs   &
         ,fac    &
         ,fac2   &
         ,dex12  &
         ,ro


    integer :: i,j,k

    ! local arrays:
    real, dimension(MAXLIMAX)::   xkh100  &
         ,xkhs    &
         ,xkdz    &
         ,xkzi    &
         ,hs


    real, dimension(MAXLIMAX,MAXLJMAX) :: help



    sm = 0.04


    xfrco=0.5*(sqrt(6859.)-1)
    exfrco=1./3.




    !c..exchange parameter and its vertical derivative at z = hs

    do j=1,ljmax
       do i=1,limax

          xkh100(i)=0.
          xkhs(i)=0.
          xkdz(i)=0.
          xkzi(i)=0.
          h100 = zs_bnd(i,j,KMAX_MID)
          !
          !
          !...................................................................
          !..air density at ground level is always calculated diagnostically:
          !

          ux0 = ustar_nwp(i,j)
          ux3 = ux0*ux0*ux0


          if(ziu(i,j) >= zimin) then
             !
             !..........................
             !..unstable surface-layer.:
             !
             !..height of surface layer
             hs(i)=sm*ziu(i,j)

             !c..hsl=hs/l where l is the monin-obhukov length
             hsl = KARMAN*GRAV*hs(i)*fh(i,j,nr)*KAPPA &
                  /(ps(i,j,nr)*ux3)


             !changes: use simple Garratt \Phi function
             !   instead of "older" Businge and Iversen/Nordeng stuff:

             xkhs(i) = ux0*KARMAN*hs(i)*sqrt(1.0-16.0*hsl)  ! /Pr=1.00
             xkdz(i) = xkhs(i)*(1.-0.5*16.0*hsl/(1.0-16.0*hsl))/hs(i)

             hsurfl = KARMAN*GRAV*h100*fh(i,j,nr)*KAPPA           &
                  /(ps(i,j,nr)*ux3)
             xkh100(i) = ux0*KARMAN*h100*sqrt(1.-16.*hsurfl)

             Kz_min(i,j)=xkh100(i)
             Kz_m2s(i,j,KMAX_MID)=xkhs(i)

     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
          write(6,"(a,4f8.2,f12.4)") "HMIX INU", ux0, ziu(i,j),  hs(i), hsl, xkhs(i)
     end if

          else
             !
             !..........................
             !..stable surface-layer...:
             !----------------------------------
             !
             !..height of surface layer
             hs(i)=sm*zixx(i,j)
             !
             !..hsl=hs/l where l is the monin-obhukov length
             hsl = KARMAN*GRAV*hs(i)*amax1(0.001,fh(i,j,nr))*KAPPA&
                  /(ps(i,j,nr)*ux3)

             Kz_m2s(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(1.00+5.0*hsl)

     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
         invL = hsl/hs(i)

!function BrostWyngaardKz(z,h,ustar,invL) result(Kz)
!and try again:
         KBW = BrostWyngaardKz( zs_bnd(i,j,KMAX_MID), zixx(i,j), ux0, invL      )
          write(6,"(a,10es12.3)") "HMIX INux0", ustar_nwp(i,j), ux0
          write(6,"(a,f9.4,3f8.2,10es12.3)") "HMIX INS", ux0, ziu(i,j), fh(i,j,nr),  hs(i), hsl, Kz_m2s(i,j,KMAX_MID), KBW
     end if


          endif

          hsurfl = KARMAN*GRAV*100.*amax1(0.001,fh(i,j,nr))*KAPPA&
               /(ps(i,j,nr)*ux3)
          Kz_min(i,j)=ux0*KARMAN*h100/(1.00+5.0*hsurfl)
     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
          write(6,"(a,4f8.2,f12.4)") "HMIX KzMIN", ux0, ziu(i,j),  hsurfl, Kz_min(i,j)
     end if
          !
          !...............................................................

       end do
       !
       !
       !..exchange parameter at z = ziu
       !
       do  k=1,KMAX_MID
          do  i=1,limax

             if(ziu(i,j).gt.zimin .and. zs_bnd(i,j,k).ge.ziu(i,j)) then
                xkzi(i)=Kz_m2s(i,j,k)
     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
         KBW = BrostWyngaardKz( zs_bnd(i,j,k), zixx(i,j), ux0, invL      )
          write(6,"(a,i3, f10.4,4f8.2,f12.4)") "HMIX KKA", k, Kz_m2s(i,j,k), ziu(i,j), zixx(i,j), invL,  KBW
     end if
             elseif (ziu(i,j).gt.zimin) then
                !
                !.....................................................
                !..the obrien-profile for z<ziu                      .
                !.....................................................
                !
                if(zs_bnd(i,j,k).le.hs(i)) then
                   Kz_m2s(i,j,k)=zs_bnd(i,j,k)*xkhs(i)/hs(i)
     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
          write(6,"(a,i3, f10.4,4f8.2,f12.4)") "HMIX KKB", k, Kz_m2s(i,j,k)
     end if
                else
                   zimhs = ziu(i,j)-hs(i)
                   zimz  =ziu(i,j)-zs_bnd(i,j,k)
                   zmhs  =zs_bnd(i,j,k)-hs(i)
                   Kz_m2s(i,j,k) = xkzi(i)+(zimz/zimhs)*(zimz/zimhs)  &
                        *(xkhs(i)-xkzi(i)+zmhs*(xkdz(i)     &
                        + 2.*(xkhs(i)-xkzi(i))/zimhs))
     if( debug_flag .and. i == debug_iloc .and. j == debug_jloc ) then
          write(6,"(a,i3, f10.4,4f8.2,f12.4)") "HMIX KKC", k, Kz_m2s(i,j,k)
     end if
                endif

             endif

          end do
       end do

    end do
    !
    !..spatial smoothing of Kz_ms2:
    !
    !

    do  k=2,KMAX_MID

       do i=1,limax
          do j=1,ljmax
             if ( (pzpbl(i,j)>z_mid(i,j,k)) )then
                Kz_m2s(i,j,k)=max(Kz_m2s(i,j,k),Kz_min(i,j))
             endif

             help(i,j) = Kz_m2s(i,j,k)
          enddo
       enddo

       call smoosp(help,KZ_MINIMUM ,KZ_MAXIMUM )

       do i=1,limax
          do j=1,ljmax
             Kz_m2s(i,j,k) = help(i,j)

             fac   = GRAV/(ps(i,j,nr) - PT)
             fac2  = fac*fac
             dex12 = th(i,j,k-1,nr)*(exnm(i,j,k) - exns(i,j,k))       &
                  + th(i,j,k,nr)*(exns(i,j,k) - exnm(i,j,k-1))
             ro    = ((ps(i,j,nr) - PT)*sigma_bnd(k) + PT)*CP*(exnm(i,j,k) &
                  - exnm(i,j,k-1))/(RGAS_KG*exns(i,j,k)*dex12)
             SigmaKz(i,j,k,nr) = Kz_m2s(i,j,k)*ro*ro*fac2
          enddo
       enddo

    end do


    !
    !...............................................................
    !..mixing-layer parameterization finished.......................
    !...............................................................
    !

  end subroutine O_Brian











  subroutine Getmeteofield(meteoname,namefield,nrec,&
       ndim,validity,field)
    !
    ! Read the meteofields and distribute to nodes
    !


    implicit none

    real, dimension(*),intent(out)  :: field ! dimensions: (MAXLIMAX,MAXLJMAX)
    ! or     (MAXLIMAX,MAXLJMAX,KMAX)

    character (len = *),intent(in)  ::meteoname,namefield
    character (len = *),intent(out) ::validity

    integer,intent(in)              :: nrec,ndim



    integer*2 :: var_local(MAXLIMAX,MAXLJMAX,KMAX_MID)
    integer*2, allocatable ::var_global(:,:,:)   ! faster if defined with
    ! fixed dimensions for all
    ! nodes?
    real :: scalefactors(2)
    integer :: KMAX,ijk,i,k,j,nfetch

    validity=''

    if(ndim==3)KMAX=KMAX_MID
    if(ndim==2)KMAX=1
    if(MasterProc)then
       allocate(var_global(GIMAX,GJMAX,KMAX))
       nfetch=1
       call GetCDF_short(namefield,meteoname,var_global,GIMAX,IRUNBEG,GJMAX, &
            JRUNBEG,KMAX,nrec,nfetch,scalefactors,validity)
    else
       allocate(var_global(1,1,1)) !just to have the array defined
    endif

    !note: var_global is defined only for me=0
    call global2local_short(var_global,var_local,MSG_READ4,GIMAX,GJMAX,&
         KMAX,1,1)
    CALL MPI_BCAST(scalefactors,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(validity,20,MPI_BYTE,0,MPI_COMM_WORLD,INFO)


    deallocate(var_global)


    ijk=0
    do k=1,KMAX ! KMAX is =1 for 2D arrays
       do j=1,MAXLJMAX
          do i=1,MAXLIMAX
             ijk=ijk+1
             field(ijk)=var_local(i,j,k)*scalefactors(1)+scalefactors(2)
          enddo
       enddo
    enddo

  end subroutine Getmeteofield







  subroutine GetCDF_short(varname,fileName,var,GIMAX,IRUNBEG,GJMAX,JRUNBEG &
       ,KMAX,nstart,nfetch,scalefactors,validity)
    !
    ! open and reads CDF file
    !
    ! The nf90 are functions which return 0 if no error occur.
    ! check is only a subroutine which check wether the function returns zero
    !
    !
    implicit none

    character (len=*),intent(in) :: fileName

    character (len = *),intent(in) ::varname
    character (len = *),intent(out) ::validity
    real,intent(out) :: scalefactors(2)
    integer, intent(in) :: nstart,GIMAX,IRUNBEG,GJMAX,JRUNBEG,KMAX
    integer, intent(inout) ::  nfetch
    integer*2, dimension(GIMAX*GJMAX*KMAX*NFETCH),intent(out) :: var
    integer :: varID,ndims
    integer :: ncFileID,var_date,status
    real :: scale,offset
    character *100 :: period_read

    validity='                                     ' !initialisation
    period_read='                                     ' !initialisation
    scalefactors(1) = 1.0 !default
    scalefactors(2) = 0.  !default

    ndims=3
    if(KMAX==1)ndims=2
    !open an existing netcdf dataset
    call check(nf90_open(path=trim(fileName),mode=nf90_nowrite,ncid=ncFileID))

    !get varID:
    status = nf90_inq_varid(ncid=ncFileID,name=trim(varname),varID=VarID)
    if(status /= nf90_noerr)then
       validity=field_not_found
       var=0
       call check(nf90_close(ncFileID))
       return
    endif


    !get scale factors

    status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
    if(status == nf90_noerr) scalefactors(1) = scale
    status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
    if(status == nf90_noerr) scalefactors(2) = offset

    !find validity
    status = nf90_get_att(ncFileID, VarID, "validity", period_read  )
    if(status == nf90_noerr)then
       validity  = trim(period_read)
    else
       status = nf90_get_att(ncFileID, VarID, "period_of_validity", &
            period_read  )
       if(status /= nf90_noerr)then
          validity='instantaneous' !default
       endif
    endif

    ! if(Nfetch<nrecords)then
    !    write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
    ! endif

    !get variable
    if(ndims==2)then
       call check(nf90_get_var(ncFileID, VarID, var,&
            start=(/IRUNBEG,JRUNBEG,nstart/),count=(/ GIMAX,GJMAX,nfetch /)))
    elseif(ndims==3)then
       call check(nf90_get_var(ncFileID, VarID, var,&
            start=(/IRUNBEG,JRUNBEG,1,nstart/),count=(/GIMAX,GJMAX,KMAX,nfetch /)))
    endif

    call check(nf90_close(ncFileID))

  end subroutine GetCDF_short





  subroutine Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm_i,xm_j,xm2,&
       ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
       ,cyclicgrid)
    !
    ! Get grid and time parameters as defined in the meteo file
    ! Do some checks on sizes and dates
    !
    ! This routine is called only once (and is therefore not optimized for speed)
    !

    implicit none

    character (len = *), intent(in) ::meteoname
    integer, intent(in):: nyear,nmonth,nday,nhour
    real, intent(out) :: GRIDWIDTH_M,xp,yp,fi, ref_latitude,&
         xm2(0:MAXLIMAX+1,0:MAXLJMAX+1)&
         ,xm_i(0:MAXLIMAX+1,0:MAXLJMAX+1)&
         ,xm_j(0:MAXLIMAX+1,0:MAXLJMAX+1),sigma_mid(KMAX_MID)
    integer, intent(out):: Nhh,nhour_first,cyclicgrid

    integer :: nseconds(1),n1,i,j,im,jm,i0,j0
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID
    integer :: GIMAX_file,GJMAX_file,KMAX_file,ihh,ndate(4)
    real,dimension(-1:GIMAX+2,-1:GJMAX+2) ::xm_global,xm_global_j,xm_global_i
    integer :: status,iglobal,jglobal,info,South_pole,North_pole,Ibuff(2)
    real :: xm_i_max,ndays(1),x1,x2,x3,x4
    character (len = 50) :: timeunit

    if(MasterProc)then
       !open an existing netcdf dataset
       status = nf90_open(path=trim(meteoname),mode=nf90_nowrite,ncid=ncFileID)

       if(status /= nf90_noerr) then
          print *,'not found',trim(meteoname)
          METEOfelt=1
       else
          print *,'  reading ',trim(meteoname)
          projection=''
         call check(nf90_get_att(ncFileID,nf90_global,"projection",projection))
          if(trim(projection)=='Rotated_Spherical'.or.trim(projection)=='rotated_spherical'&
               .or.trim(projection)=='rotated_pole'.or.trim(projection)=='rotated_latitude_longitude')then
             projection='Rotated_Spherical'
          endif
          write(*,*)'projection: ',trim(projection)

          !get dimensions id
          if(trim(projection)=='Stereographic') then
             call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
          elseif(trim(projection)==trim('lon lat')) then
             call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
          else
             !     write(*,*)'GENERAL PROJECTION ',trim(projection)
             call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
          endif

          call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))
          call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))

          !get dimensions length
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))

          write(*,*)'dimensions meteo grid',GIMAX_file,GJMAX_file,KMAX_file,Nhh

          if(GIMAX_file/=IIFULLDOM.or.GJMAX_file/=JJFULLDOM)then
             write(*,*)'IIFULLDOM,JJFULLDOM ',IIFULLDOM,JJFULLDOM
             write(*,*)'WARNING: large domain and meteorology file have different sizes'
             write(*,*)'WARNING: THIS CASE IS NOT TESTED. Please change large domain'
          endif

          call CheckStop(GIMAX+IRUNBEG-1 > GIMAX_file, "NetCDF_ml: I outside domain" )
          call CheckStop(GJMAX+JRUNBEG-1 > GJMAX_file, "NetCDF_ml: J outside domain" )
          call CheckStop(KMAX_MID > KMAX_file,  "NetCDF_ml: wrong vertical dimension")
          call CheckStop(24/Nhh, METSTEP,          "NetCDF_ml: METSTEP != meteostep" )

          call CheckStop(nhour/=0 .and. nhour /=3,&
               "Met_ml/GetCDF: must start at nhour=0 or 3")

          call check(nf90_get_att(ncFileID,timeVarID,"units",timeunit))

          ihh=1
          n1=1
          if(trim(timeunit)==trim("days since 1900-1-1 0:0:0"))then
             write(*,*)'Meteo date in days since 1900-1-1 0:0:0'
             call check(nf90_get_var(ncFileID,timeVarID,ndays,&
                  start=(/ihh/),count=(/n1 /)))

             call datefromdayssince1900(ndate,ndays(1),0)
          else
             call check(nf90_get_var(ncFileID,timeVarID,nseconds,&
                  start=(/ihh/),count=(/n1 /)))
             call datefromsecondssince1970(ndate,nseconds(1),0) !default
          endif
          nhour_first=ndate(4)


          call CheckStop(ndate(1), nyear,  "NetCDF_ml: wrong meteo year" )
          call CheckStop(ndate(2), nmonth, "NetCDF_ml: wrong meteo month" )
          call CheckStop(ndate(3), nday,   "NetCDF_ml: wrong meteo day" )

          do ihh=1,Nhh

             if(trim(timeunit)==trim("days since 1900-1-1 0:0:0"))then
                call check(nf90_get_var(ncFileID, timeVarID, ndays,&
                     start=(/ ihh /),count=(/ n1 /)))
                call datefromdayssince1900(ndate,ndays(1),0)
             else
                call check(nf90_get_var(ncFileID, timeVarID, nseconds,&
                  start=(/ ihh /),count=(/ n1 /)))
                call datefromsecondssince1970(ndate,nseconds(1),0)
             endif

             call CheckStop( mod((ihh-1)*METSTEP+nhour_first,24), ndate(4),  &
                  "NetCDF_ml: wrong meteo hour" )

          enddo


          !get global attributes
          call check(nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M))
          if(trim(projection)=='Stereographic')then
             call check(nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude))
             call check(nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole" &
                  ,xp ))
             call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole" &
                  ,yp ))
             call check(nf90_get_att(ncFileID, nf90_global, "fi",fi ))

             call GlobalPosition
          elseif(trim(projection)==trim('lon lat')) then
             ref_latitude=60.
             xp=0.0
             yp=GJMAX
             fi =0.0
             call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
             call check(nf90_get_var(ncFileID, varID, gl_glob(1:IIFULLDOM,1) ))
             do i=1,IIFULLDOM
                if(gl_glob(i,1)>180.0)gl_glob(i,1)=gl_glob(i,1)-360.0
                if(gl_glob(i,j)<-180.0)gl_glob(i,j)=gl_glob(i,j)+360.0
             enddo
             do j=1,JJFULLDOM
                gl_glob(:,j)=gl_glob(:,1)
             enddo
             call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
             call check(nf90_get_var(ncFileID, varID, gb_glob(1,1:JJFULLDOM) ))
             do i=1,IIFULLDOM
                gb_glob(i,:)=gb_glob(1,:)
             enddo
          else
             ref_latitude=60.
             xp=0.0
             yp=GJMAX
             fi =0.0
             if(trim(projection)=='Rotated_Spherical')then
                call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_latitude",grid_north_pole_latitude))
                call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_longitude",grid_north_pole_longitude))
             endif
             call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
             call check(nf90_get_var(ncFileID, varID, gl_glob(1:IIFULLDOM,1:JJFULLDOM) ))

             call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
             call check(nf90_get_var(ncFileID, varID, gb_glob(1:IIFULLDOM,1:JJFULLDOM) ))
             do j=1,JJFULLDOM
             do i=1,IIFULLDOM
                if(gl_glob(i,j)>180.0)gl_glob(i,j)=gl_glob(i,j)-360.0
                if(gl_glob(i,j)<-180.0)gl_glob(i,j)=gl_glob(i,j)+360.0
             enddo
             enddo

          endif
          !get variables
          status=nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID)

          if(status == nf90_noerr)then
             !mapping factor at center of cells is defined
             !make "staggered" map factors
             call check(nf90_get_var(ncFileID, varID, xm_global(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
             do j=1,GJMAX
                do i=1,GIMAX-1
                   xm_global_j(i,j)=0.5*(xm_global(i,j)+xm_global(i+1,j))
                enddo
             enddo
             i=GIMAX
             do j=1,GJMAX
                xm_global_j(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i-1,j)
             enddo
             do j=1,GJMAX-1
                do i=1,GIMAX
                   xm_global_i(i,j)=0.5*(xm_global(i,j)+xm_global(i,j+1))
                enddo
             enddo
             j=GJMAX
             do i=1,GIMAX
                xm_global_i(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i,j-1)
             enddo

          else
             !map factor are already staggered
             status=nf90_inq_varid(ncid=ncFileID, name="map_factor_i", varID=varID)

             !BUGCHECK - moved here... (deleted if loop)
             call CheckStop( status, nf90_noerr, "erro rreading map factor" )

             call check(nf90_get_var(ncFileID, varID, xm_global_i(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
             call check(nf90_inq_varid(ncid=ncFileID, name="map_factor_j", varID=varID))
             call check(nf90_get_var(ncFileID, varID, xm_global_j(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
          endif

          call check(nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID))
          call check(nf90_get_var(ncFileID, varID, sigma_mid ))


       endif !found meteo
    endif !me=0

    CALL MPI_BCAST(METEOfelt ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    if( METEOfelt==1)return !do not use NetCDF meteo input


    CALL MPI_BCAST(Nhh,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(GRIDWIDTH_M,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(ref_latitude,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(yp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(fi,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(sigma_mid,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_i(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_j(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(gb_glob(1:IIFULLDOM,1:JJFULLDOM),8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(gl_glob(1:IIFULLDOM,1:JJFULLDOM),8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(projection,len(projection),MPI_CHARACTER,0,MPI_COMM_WORLD,INFO) 


    do j=1,MAXLJMAX
       do i=1,MAXLIMAX
          gl(i,j)=gl_glob(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          gb(i,j)=gb_glob(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
       enddo
    enddo
    i0=0
    im=MAXLIMAX
    j0=0
    jm=MAXLJMAX
    if(gi0+MAXLIMAX+1+IRUNBEG-2>IIFULLDOM)im=MAXLIMAX-1!outside fulldomain
    if(gi0+0+IRUNBEG-2<1)i0=1!outside fulldomain
    if(gj0+MAXLJMAX+1+JRUNBEG-2>JJFULLDOM)jm=MAXLJMAX-1!outside fulldomain
    if(gj0+0+JRUNBEG-2<1)j0=1!outside fulldomain

    do j=j0,jm
       do i=i0,im
          x1=gl_glob(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x2=gl_glob(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x3=gl_glob(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)
          x4=gl_glob(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2)

!8100=90*90; could use any number much larger than zero and much smaller than 180*180  
          if(x1*x2<-8100.0 .or. x1*x3<-8100.0 .or. x1*x4<-8100.0)then
             !Points are on both sides of the longitude -180=180              
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
             if(x3<0)x3=x3+360.0
             if(x4<0)x4=x4+360.0
          endif
          gl_stagg(i,j)=0.25*(x1+x2+x3+x4)

          gb_stagg(i,j)=0.25*(gb_glob(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               gb_glob(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               gb_glob(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
               gb_glob(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
       enddo
    enddo
    do j=0,j0
       do i=i0,im
          x1=gl_stagg(i,j+1)
          x2=gl_stagg(i,j+2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j+1)-gb_stagg(i,j+2)
       enddo
    enddo
    do j=jm,MAXLJMAX
       do i=i0,im
          x1=gl_stagg(i,j-1)
          x2=gl_stagg(i,j-2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j-1)-gb_stagg(i,j-2)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=0,i0
          x1=gl_stagg(i+1,j)
          x2=gl_stagg(i+2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i+1,j)-gb_stagg(i+2,j)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=im,MAXLIMAX
          x1=gl_stagg(i-1,j)
          x2=gl_stagg(i-2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i-1,j)-gb_stagg(i-2,j)
       enddo
    enddo
!ensure that values are within [-180,+180]]
    do j=0,MAXLJMAX
       do i=0,MAXLIMAX
          if(gl_stagg(i,j)>180.0)gl_stagg(i,j)=gl_stagg(i,j)-360.0
          if(gl_stagg(i,j)<-180.0)gl_stagg(i,j)=gl_stagg(i,j)+360.0
       enddo
    enddo

    !test if the grid is cyclicgrid:
    !The last cell + 1 cell = first cell
    Cyclicgrid=1 !Cyclicgrid
    do j=1,JJFULLDOM
       if(mod(nint(gl_glob(GIMAX,j)+360+360.0/GIMAX),360)/=&
            mod(nint(gl_glob(IRUNBEG,j)+360.0),360))then
          Cyclicgrid=0  !not cyclicgrid
       endif
    enddo

    if(MasterProc .and. DEBUG_MET)write(*,*)'CYCLICGRID:',Cyclicgrid

    !complete (extrapolate) along the four lateral sides
    do i=1,GIMAX
       xm_global_j(i,0)=1.0/(2.0/(xm_global_j(i,1))-1.0/(xm_global_j(i,2)))
       xm_global_j(i,-1)=1.0/(2.0/(xm_global_j(i,0))-1.0/(xm_global_j(i,1)))
       xm_global_j(i,GJMAX+1)=1.0/(2.0/(xm_global_j(i,GJMAX))-1.0/(xm_global_j(i,GJMAX-1)))
       xm_global_j(i,GJMAX+2)=1.0/(2.0/(xm_global_j(i,GJMAX+1))-1.0/(xm_global_j(i,GJMAX)))
       xm_global_i(i,0)=1.0/(2.0/(xm_global_i(i,1))-1.0/(xm_global_i(i,2)))
       xm_global_i(i,-1)=1.0/(2.0/(xm_global_i(i,0))-1.0/(xm_global_i(i,1)))
       xm_global_i(i,GJMAX+1)=1.0/(2.0/(xm_global_i(i,GJMAX))-1.0/(xm_global_i(i,GJMAX-1)))
       xm_global_i(i,GJMAX+2)=1.0/(2.0/(xm_global_i(i,GJMAX+1))-1.0/(xm_global_i(i,GJMAX)))
    enddo
    do j=-1,GJMAX+2
       xm_global_j(0,j)=1.0/(2.0/(xm_global_j(1,j))-1.0/(xm_global_j(2,j)))
       xm_global_j(-1,j)=1.0/(2.0/(xm_global_j(0,j))-1.0/(xm_global_j(1,j)))
       xm_global_j(GIMAX+1,j)=1.0/(2.0/(xm_global_j(GIMAX,j))-1.0/(xm_global_j(GIMAX-1,j)))
       xm_global_j(GIMAX+2,j)=1.0/(2.0/(xm_global_j(GIMAX+1,j))-1.0/(xm_global_j(GIMAX,j)))
       xm_global_i(0,j)=1.0/(2.0/(xm_global_i(1,j))-1.0/(xm_global_i(2,j)))
       xm_global_i(-1,j)=1.0/(2.0/(xm_global_i(0,j))-1.0/(xm_global_i(1,j)))
       xm_global_i(GIMAX+1,j)=1.0/(2.0/(xm_global_i(GIMAX,j))-1.0/(xm_global_i(GIMAX-1,j)))
       xm_global_i(GIMAX+2,j)=1.0/(1.0/(2*xm_global_i(GIMAX+1,j))-1.0/(xm_global_i(GIMAX,j)))
    enddo

    j=1
    i=1
    if(abs(1.5*gb_glob(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)-0.5*gb_glob(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2))>89.5)then
       write(*,*)'south pole' !xm is infinity
       xm_global_i(:,0)=1.0E19
       xm_global_i(:,-1)=1.0E19
    endif



    !keep only part of xm relevant to the local domain
    !note that xm has dimensions larger than local domain

    call CheckStop( MAXLIMAX+1 > limax+2, "Error in Met_ml X size definition" )
    call CheckStop( MAXLJMAX+1 > ljmax+2, "Error in Met_ml J size definition" )

    do j=0,MAXLJMAX+1
       do i=0,MAXLIMAX+1
          iglobal=gi0+i-1
          jglobal=gj0+j-1
          xm_i(i,j)=xm_global_i(iglobal,jglobal)
          xm_j(i,j)=xm_global_j(iglobal,jglobal)
          !Note that xm is inverse length: interpolate 1/xm rather than xm
          xm2(i,j) = 4.0*( (xm_global_i(iglobal,jglobal-1)*&
               xm_global_i(iglobal,jglobal))/ &
               (xm_global_i(iglobal,jglobal-1)+&
               xm_global_i(iglobal,jglobal)    )   )   *(&
               xm_global_j(iglobal-1,jglobal)*&
               xm_global_j(iglobal,jglobal)     )/(&
               xm_global_j(iglobal-1,jglobal)+&
               xm_global_j(iglobal,jglobal)     )
       enddo
    enddo

    !pw
    !If some cells are to narrow (Poles in lat lon coordinates),
    !this will give too small time steps in the Advection,
    !because of the constraint that the Courant number should be <1.
    !
    !If Poles are found and lon-lat coordinates are used the Advection scheme
    !will be modified to be able to cope with the singularity

    !Look for poles
    !If the northernmost or southernmost lines are poles, they are not
    !considered as outer boundaries and will not be treated 
    !by "BoundaryConditions_ml".
    !Note that "Poles" is defined in subdomains

    North_pole=1
    do i=1,limax
       if(nint(gb(i,ljmax))<=88)then
          North_pole=0  !not north pole
       endif
    enddo

    South_pole=1
    do i=1,limax
       if(nint(gb(i,1))>=-88)then
          South_pole=0  !not south pole
       endif
    enddo

    Poles=0
    if(North_pole==1)then
       Poles(1)=1
       write(*,*)me,'Found North Pole'
    endif

    if(South_pole==1)then
       Poles(2)=1
       write(*,*)me,'Found South Pole'
    endif

     CALL MPI_ALLREDUCE(Poles,Ibuff,2,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,INFO)
     Pole_included=max(Ibuff(1),Ibuff(2))
     
     if(ME==0.and.Pole_included==1)write(*,*)'The grid includes a pole'
     
  end subroutine Getgridparams





  subroutine check(status)
    implicit none
    integer, intent ( in) :: status

    call CheckStop( status, nf90_noerr, "Error in Met_ml/NetCDF stuff" &
         //  trim( nf90_strerror(status) ) )

  end subroutine check
  !_______________________________________________________________________

  subroutine datefromsecondssince1970(ndate,nseconds,printdate)
    ! calculate date from seconds that have passed since the start of
    ! the year 1970

    implicit none

    integer, intent(out) :: ndate(4)
    integer, intent(in) :: nseconds
    integer,  intent(in) :: printdate

    integer :: n,nday,nmdays(12),nmdays2(13)
    nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    nmdays2(1:12)=nmdays
    nmdays2(13)=0
    ndate(1)=1969
    n=0
    do while(n<=nseconds)
       n=n+24*3600*365
       ndate(1)=ndate(1)+1
       if(mod(ndate(1),4)==0)n=n+24*3600
    enddo
    n=n-24*3600*365
    if(mod(ndate(1),4)==0)n=n-24*3600
    if(mod(ndate(1),4)==0)nmdays2(2)=29
    ndate(2)=0
    do while(n<=nseconds)
       ndate(2)=ndate(2)+1
       n=n+24*3600*nmdays2(ndate(2))
    enddo
    n=n-24*3600*nmdays2(ndate(2))
    ndate(3)=0
    do while(n<=nseconds)
       ndate(3)=ndate(3)+1
       n=n+24*3600
    enddo
    n=n-24*3600
    ndate(4)=-1
    do while(n<=nseconds)
       ndate(4)=ndate(4)+1
       n=n+3600
    enddo
    n=n-3600
    !    ndate(5)=nseconds-n
    if(printdate>0)then
       write(*,*)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
            ndate(3),', hour: ',ndate(4),', seconds: ',nseconds-n
    endif
  end subroutine datefromsecondssince1970


  subroutine datefromdayssince1900(ndate,ndays,printdate)
    ! calculate date from days that have passed since the start of
    ! the year 1900
    ! NB: 1900 is not a leap year


    implicit none

    integer, intent(out) :: ndate(4)
    real, intent(in) :: ndays
    integer,  intent(in) :: printdate

    integer :: n,nday,nmdays(12),nmdays2(13)
    real*8 :: nr

    nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    nmdays2(1:12)=nmdays
    nmdays2(13)=0
    ndate(1)=1899
    n=0
    do while(n<=ndays)
       n=n+365
       ndate(1)=ndate(1)+1
       if(mod(ndate(1),4)==0.and.ndate(1)/=1900)n=n+1 !NB: 1900 is not a leap year
    enddo
    n=n-365
    if(mod(ndate(1),4)==0.and.ndate(1)/=1900)n=n-1
    if(mod(ndate(1),4)==0.and.ndate(1)/=1900)nmdays2(2)=29
    ndate(2)=0
    do while(n<=ndays)
       ndate(2)=ndate(2)+1
       n=n+nmdays2(ndate(2))
    enddo
    n=n-nmdays2(ndate(2))
    ndate(3)=0
    do while(n<=ndays)
       ndate(3)=ndate(3)+1
       n=n+1
    enddo
    n=n-1
    ndate(4)=-1
    nr=n
    do while(nr<=ndays+1.0d0/24.0d0/3600.0d0/2)
       ndate(4)=ndate(4)+1
       nr=nr+1.0d0/24.0d0
    enddo
    nr=nr-1.0d0/24.0d0
    if(nint((ndays-nr)*3600.0*24.0)/=0)then
       write(*,*)'WARNING: not integer number of hours'
    endif
    !    ndate(5)=nseconds-n
    if(printdate>0)then
       write(*,*)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
            ndate(3),', hour: ',ndate(4),', seconds: ',nint((ndays-nr)*3600.0*24.0)
    endif
  end subroutine datefromdayssince1900

end module met_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



