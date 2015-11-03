subroutine hourly_out() !!  spec,ofmt,ix1,ix2,iy1,iy2,unitfac)
!***********************************************************************
!**    DESCRIPTION:
!       Calculates and
!       Outputs hourly concentration (or met) values for a sub-set of the grid.
!
!** Surface output
! Onyl relevant for lowermost model level, meaningless for higher levels.
! - ADVppbv: instantaneous surface concentrations in ppb.
! - ADVugXX: instantaneous surface concentrations in ug (ug/m3, ugC/m3, ugS/m3, ugN/m3).
!   For ug/m3     output, set hr_out%unitconv=to_ug_ADV(ixadv).
!   For ugX/m3    output, set hr_out%unitconv=to_ug_X(ixadv).
! - D2D_mean: hourly means comparable to daily output defined on My_Derived_ml.
! - D2D_inst: instantaneous value used in daily output defined on My_Derived_ml.
! - D2D_accum: accumulated   alue used in daily output defined on My_Derived_ml.
! - D2D: D2D_mean/D2D_accum according to the corresponding My_Derived_ml definition.
!
!** Column integrated output
! - COLUMN: 
!   For ug/m2     output, set hr_out%unitconv=to_ug_ADV(ixadv).
!   For ugX/m2    output, set hr_out%unitconv=to_ug_X(ixadv).
!   For molec/cm2 output, set hr_out%unitconv=to_molec_cm2.
!
!** Multi-layer output
!  NLEVELS_HOURLY (My_Outputs_ml) max levels
!  hr_out%nk levels for each output
! - BCVppbv: instantaneous grid-centre concentrations in ppb.
! - BCVugXX: instantaneous grid-centre concentrations in ug (ug/m3, ugC/m3, ugS/m3, ugN/m3).
!   For ug/m3     output, set hr_out%unitconv=to_ug_ADV(ixadv).
!   For ugX/m3    output, set hr_out%unitconv=to_ug_X(ixadv).
!*************************************************************************
use My_Outputs_ml,    only: NHOURLY_OUT,    & ! No. outputs
                            NLEVELS_HOURLY, & ! No. output levels
                            hr_out,         & ! Required outputs
                            LEVELS_HOURLY ! Output selected model levels
          !NML SELECT_LEVELS_HOURLY, LEVELS_HOURLY ! Output selected model levels

use CheckStop_ml,     only: CheckStop
use Chemfields_ml,    only: xn_adv,xn_shl,cfac,PM25_water,PM25_water_rh50
use ChemGroups_ml,    only: chemgroups
use Derived_ml,       only: num_deriv2d,nav_2d        ! D2D houtly output type
use DerivedFields_ml, only: f_2d,d_2d          ! D2D houtly output type
use OwnDataTypes_ml,  only: Asc2D, Deriv
use ChemSpecs,        only: NSPEC_SHL, species
!CMR  use ChemSpecs_shl_ml ,only: NSPEC_SHL          ! Maps indices
!CMR  use ChemChemicals_ml ,only: species            ! Gives names
use GridValues_ml,    only: i_fdom, j_fdom,&   ! Gives emep coordinates
                            debug_proc, debug_li,debug_lj
use Io_ml,            only: IO_HOURLY
use ModelConstants_ml,only: KMAX_MID, MasterProc, MY_OUTPUTS,&
                            IOU_INST, IOU_HOUR, IOU_YEAR, IOU_YEAR_LASTHH, &
                            DEBUG => DEBUG_OUT_HOUR,runlabel1,HOURLYFILE_ending,&
                            FORECAST
use ModelConstants_ml,only: SELECT_LEVELS_HOURLY !NML
use MetFields_ml,     only: t2_nwp,th, roa, surface_precip, ws_10m ,rh2m,&
                            pzpbl, ustar_nwp, Kz_m2s, &
                            Idirect, Idiffuse, z_bnd, z_mid,ps
use NetCDF_ml,        only: Out_netCDF, CloseNetCDF, Init_new_netCDF, fileName_hour, &
                            Int1, Int2, Int4, Real4, Real8  !Output data type to choose
use OwnDataTypes_ml,  only: TXTLEN_DERIV,TXTLEN_SHORT
use Par_ml,           only: LIMAX, LJMAX, GIMAX,GJMAX,        &
                            me, IRUNBEG, JRUNBEG, limax, ljmax
use Pollen_ml,        only: heatsum, pollen_left, AreaPOLL
use TimeDate_ml,      only: current_date
use TimeDate_ExtraUtil_ml,only : date2string
use Units_ml,         only: Group_Units,&
                            to_number_cm3 ! converts roa [kg/m3] to M [molec/cm3]

implicit none

!*.. Components of  hr_out
!* character(len=TXTLEN_DERIV):: name   ! netCDF variable name
!* character(len=TXTLEN_SHORT):: type   ! "ADVppbv" or "ADVugm3" or "SHLmcm3"
!* character(len=9) :: ofmt     ! Output format (e.g. es12.4)
!* integer          :: spec     ! Species number in xn_adv or xn_shl array or other arrays
!* integer          :: ix1,ix2  ! bottom-left & upper-right x
!* integer          :: iy1,iy2  ! bottom-left & upper-right y
!* integer          :: nk       ! number of vertical levels
!* character(len=TXTLEN_SHORT) :: unit   ! netCDF unit attribute
!* real             :: unitconv !  conv. factor
!* real             :: max      ! Max allowed value for output

! local variables
  logical, save     :: first_call = .true. ! Set false after file opened
  real hourly(LIMAX,LJMAX)      ! Local hourly value  (e.g. ppb)
  real ghourly(GIMAX,GJMAX)           ! Global hourly value (e.g. ppb)
  real :: arrmax                      ! Maximum value from array
  real :: unit_conv                   ! Unit conversion (ppb ug etc.)
  real :: g                           ! tmp - saves value of ghourly(i,j)
  integer, dimension(2) :: maxpos     ! Location of max value
  integer i,j,ih,ispec,itot,iadv      ! indices
  integer :: k,ik,iik                 ! Index for vertical level
  integer ist,ien,jst,jen             ! start and end coords
  character(len=TXTLEN_DERIV) :: name ! For output file, species names
  character(len=4)  :: suffix         ! For date "mmyy"
  integer, save :: prev_month = -99   ! Initialise with non-possible month
  type(Deriv) :: def1           ! for NetCDF
  real :: scale                 ! for NetCDF
  integer ::CDFtype,nk,klevel,ncfileID   ! for NetCDF
  character(len=TXTLEN_SHORT)    :: hr_out_type=""      ! hr_out%type
  integer                        :: hr_out_nk=0         ! hr_out%nk
  integer, pointer, dimension(:) :: gspec=>null()       ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null()  ! & unit conv. factors

  integer, allocatable, dimension(:), save :: navg ! D2D average counter

  character(len=len(fileName_hour))      :: filename
  character(len=len(fileName_hour)),save :: filename_old="no file"
  logical, save :: debug_flag      ! = ( MasterProc .and. DEBUG )
  logical       :: surf_corrected=.true.  ! to get 3m values

  character(len=*), parameter :: &
    SRF_TYPE(9)=(/"ADVppbv     ","ADVugXX     ","ADVugXXgroup",&
                  "COLUMN      ","COLUMNgroup ","D2D         ",&
                  "D2D_mean    ","D2D_inst    ","D2D_accum   "/)

  character(len=52) :: errmsg = "ok"

  if(NHOURLY_OUT<=0)then
    if(first_call.and.MasterProc.and.DEBUG) &
      write(*,*)"DEBUG Hourly_out: nothing to output!"
    first_call=.false.
    return
  endif

  ! only write at 12UTC for "TRENDS@12UTC", eg
  ! if(MY_OUTPUTS=="TRENDS@12UTC".and.current_date%hour/=12)return
  i=index(MY_OUTPUTS,"@")
  if(i>0)then
    read(MY_OUTPUTS(i+1:i+2),*)ih
    if(current_date%hour/=ih)return
  endif

  if(first_call) then
    first_call = .false.

    debug_flag=(debug_proc.and.DEBUG)
  !/ Ensure that domain limits specified in My_Outputs lie within
  !  model domain. In emep coordinates we have:
    do ih = 1, NHOURLY_OUT
      hr_out(ih)%ix1 = max(IRUNBEG,hr_out(ih)%ix1)
      hr_out(ih)%iy1 = max(JRUNBEG,hr_out(ih)%iy1)
      hr_out(ih)%ix2 = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
      hr_out(ih)%iy2 = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
!     hr_out(ih)%nk  = min(KMAX_MID,hr_out(ih)%nk)
      if(debug_flag) write(*,*) "DEBUG Hourly nk ", ih, hr_out(ih)%nk
    enddo ! ih

    allocate(navg(NHOURLY_OUT)) ! allocate and initialize
    navg(:)=0.0                 ! D2D average counter
  endif  ! first_call

  filename=trim(runlabel1)//date2string(trim(HOURLYFILE_ending),current_date)
  if(filename/=filename_old)then
    filename_old=filename

    if(debug_flag)&
      write(*,*) "DEBUG ",trim(HOURLYFILE_ending),"-Hourlyfile ",trim(filename)
  !! filename will be overwritten
    call Init_new_netCDF(trim(filename),IOU_HOUR)

    ncfileID=-1 ! must be <0 as initial value

  !! Create variables first, without writing them (for performance purposes)   
    do ih=1,NHOURLY_OUT
      def1%name=hr_out(ih)%name
      def1%unit=hr_out(ih)%unit
      def1%class=hr_out(ih)%type
      ist = hr_out(ih)%ix1
      jst = hr_out(ih)%iy1
      ien = hr_out(ih)%ix2
      jen = hr_out(ih)%iy2
      nk  = hr_out(ih)%nk
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.
      if(any(hr_out(ih)%type==SRF_TYPE))nk=1
      select case(nk)
      case(1)       ! write as 2D
        call Out_netCDF(IOU_HOUR,def1,2,1,hourly,scale,CDFtype,ist,jst,ien,jen,&
          create_var_only=.true.,ncFileID_given=ncFileID)
      case(2:)      ! write as 3D
        call Out_netCDF(IOU_HOUR,def1,3,1,hourly,scale,CDFtype,ist,jst,ien,jen,1,&
          create_var_only=.true.,chunksizes=(/ien-ist+1,jen-jst+1,1,1/),ncFileID_given=ncFileID)
      endselect
    enddo
  endif
!......... Uses concentration/met arrays from Chem_ml or Met_ml ..................
!
!        real xn_adv(NSPEC_ADV,LIMAX,LJMAX,KMAX_MID)
!        real cfac(NSPEC_ADV,LIMAX,LJMAX)
! or...
!        real xn_shl(NSPEC_ADV,LIMAX,LJMAX,KMAX_MID)
! or...
!        real temp2m(LIMAX,LJMAX)
!
!..........................................................................

  hourly(:,:) = 0.0 ! initialize

  HLOOP: do ih = 1, NHOURLY_OUT
    hr_out_nk=hr_out(ih)%nk
    if(any(hr_out(ih)%type==SRF_TYPE))hr_out_nk=1

    KVLOOP: do k = 1,hr_out_nk
      ispec = hr_out(ih)%spec
      name  = hr_out(ih)%name
      hr_out_type=hr_out(ih)%type
      if(debug_flag) &
        write(*,'(a,2i4,1X,a,/a,1X,2a,i3)')"DEBUG DERIV HOURLY", ih, ispec, &
          trim(name),"INTO HOUR TYPE:", &
          trim(hr_out(ih)%type) // " "//trim(hr_out(ih)%name), " nk:", hr_out_nk

      select case (hr_out_type)
      case("COLUMN","COLUMNgroup")
        ik=KMAX_MID-hr_out(ih)%nk+1  ! top of the column
        if(ik>=KMAX_MID)ik=1         ! 1-level column does not make sense
      case("D2D")
        ik=KMAX_MID                  ! surface/lowermost level
        if(f_2d(ispec)%avg)then ! averaged variables        
!         hr_out_type="D2D_inst"   ! output instantaneous values
          hr_out_type="D2D_mean"   ! output mean values
        else                    ! accumulated variables
          hr_out_type="D2D_accum"
        endif
      case("ADVppbv","ADVugXX","ADVugXXgroup","PMwaterSRF",&
           "D2D_inst","D2D_mean","D2D_accum")
        ik=KMAX_MID                  ! surface/lowermost level
      case default
        ik=KMAX_MID-k+1              ! all levels from model bottom are outputed,
        if(debug_flag) write(*,*)"SELECT LEVELS? ", ik, SELECT_LEVELS_HOURLY
        if(SELECT_LEVELS_HOURLY)then ! or the output levels are taken
          ik=LEVELS_HOURLY(k)        ! from LEVELS_HOURLY array (default)
          if(debug_flag) write(*,*)"DEBUG SELECT LEVELS", ik, hr_out_type
          surf_corrected = (ik==0)   ! Will implement cfac
          if(debug_flag.and.surf_corrected) &
            write(*,*)"DEBUG HOURLY Surf_correction", ik, k
!TESTHH QUERY: see below
          select case(ik)
          case(1:)
            ik=KMAX_MID-ik+1         ! model level to be outputed
          case(0)
            ik=KMAX_MID              ! surface/lowermost level
            if(debug_flag) write(*,*)"DEBUG LOWEST LEVELS", ik, hr_out_type
            select case(hr_out_type) ! ensure surface output
            case("BCVppbv","BCVugXX","BCVugXXgroup")
              hr_out_type(1:3)="ADV"
            case("PMwater")
              hr_out_type=trim(hr_out_type)//"SRF"
            endselect
          endselect
        endif
      endselect


      if(debug_flag) write(*,"(5a,i4)") "DEBUG Hourly MULTI ",&
              trim(hr_out(ih)%name), " case ", trim(hr_out_type), " k: ",  ik
      OPTIONS: select case(hr_out_type)

      case("ADVppbv")
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                      * cfac(ispec,i,j) &    ! 50m->3m conversion
                      * unit_conv            ! Units conv.
        endforall

      case("Out3D")
        itot = ispec
        iadv  = ispec - NSPEC_SHL
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        if(itot < NSPEC_SHL .and.  debug_proc) write(*,"(a,a,4es12.3)") &
          "OUT3D MAR22 ", trim(name), unit_conv, &
           to_number_cm3,roa(2,2,ik,1)*to_number_cm3, xn_shl(ispec,2,2,ik)

        !MAR22 Added SHL option.
        if( itot <= NSPEC_SHL ) then !TESTSHL
         ! CRUDE units fix. Sort out later.
         !   Inverse of No. air mols/cm3 = 1/M
         ! where M =  roa (kgair m-3) * to_number_cm3  when ! scale in ug,  else 1
         !inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * to_number_cm3 )

          forall(i=1:limax,j=1:ljmax)
             hourly(i,j) = xn_shl(ispec,i,j,ik)
          end forall
          if(index(hr_out(ih)%unit,"ppt")>0) then
            forall(i=1:limax,j=1:ljmax)
                hourly(i,j) = hourly(i,j)/( roa(i,j,ik,1) * to_number_cm3 ) * 1.0e9
            end forall
          else
            call CheckStop("SHL Out3D option not coded yet")
          end if
        
        else  ! ADV:, original code
  
          if(index(hr_out(ih)%unit,"ppb")>0) then
            forall(i=1:limax,j=1:ljmax)
              hourly(i,j) = xn_adv(iadv ,i,j,ik) &
                          * unit_conv            ! Units conv.
            end forall
          else if(index(hr_out(ih)%unit,"ug")>0) then
            forall(i=1:limax,j=1:ljmax)
              hourly(i,j) = xn_adv(iadv ,i,j,ik) &
                        * roa(i,j,ik,1)        & ! density.
                        * unit_conv              ! Units conv.
            end forall
          else
            call CheckStop("ERROR: Output_hourly  unit problem"//trim(name) )
          endif
        end if ! MAR22 ADV/SHL split 

        if(surf_corrected.and.ik==KMAX_MID.and.itot>NSPEC_SHL) then
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)*cfac(iadv,i,j) ! 50m->3m conversion
          endforall
        endif

        if(debug_flag) then
          i=debug_li; j=debug_lj
          write(*,'(A,2I4,1X,L2,2f10.4)')"Out3D K-level"//trim(name), ik,  &
            itot, surf_corrected, hourly(i,j), cfac(ispec-NSPEC_SHL,i,j)
        endif

      case("BCVppbv")
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall(i=1:limax,j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,ik) &
                      * unit_conv            ! Units conv.
        endforall
        if(DEBUG) &
          write(*,'(A,I0,1X,L2)')"K-level", ik, trim(name), itot, surf_corrected

      case("ADVugXX")  !ug/m3, ugX/m3 output at the surface
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                      * cfac(ispec,i,j) &     ! 50m->3m conversion
                      * unit_conv       &     ! Units conv.
                      * roa(i,j,KMAX_MID,1)   ! density.
        endforall

      case("BCVugXX")  ! ug/m3, ugX/m3 output at model mid-levels
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,ik) &
                !BCV  * cfac(ispec,i,j) &     ! 50m->3m conversion
                      * unit_conv       &     ! Units conv.
                      * roa(i,j,ik,1)         ! density.
        endforall
        if(DEBUG) &
          write(*,'(a,i5,a10,i5)')"K-level", ik, trim(name), itot

      case("ADVugXXgroup")  ! GROUP output in ug/m3, ugX/m3 at the surface
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = dot_product(xn_adv(gspec,i,j,KMAX_MID),&
                                    cfac(gspec,i,j) & ! 50m->3m conversion
                                   *gunit_conv(:))  & ! Units conv.
                      * roa(i,j,KMAX_MID,1)           ! density.
        endforall
        if(DEBUG) &
          write(*,'(2a10,99i5)')"Surface", trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("BCVugXXgroup")  ! GROUP output in ug/m3, ugX/m3 at model mid-levels
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = dot_product(xn_adv(gspec,i,j,ik), &
                                    gunit_conv(:))  & ! Units conv.
                      * roa(i,j,ik,1)                 ! density.
        endforall
        if(DEBUG) &
          write(*,'(a10,i7,a10,i7)')"K-level", ik, trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("PMwater")    ! PM water content in ug/m3 at model mid-levels
        if(hr_out(ih)%unit/="ug/m3")hr_out(ih)%unit="ug"
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = PM25_water(i,j,ik)

      case("PMwaterSRF")  ! PM water content in ug/m3 at surface level
        if(hr_out(ih)%unit/="ug/m3")hr_out(ih)%unit="ug"
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = PM25_water_rh50(i,j)

      case("Z","Z_MID")
        name = "Z_MID"
        unit_conv =  hr_out(ih)%unitconv
        if(surf_corrected)then
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = 3.0*unit_conv
        else
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = z_mid(i,j,ik)*unit_conv
        endif

      case("dZ","dZ_BND")  ! level thickness
        name = "dZ_BND"
        unit_conv =  hr_out(ih)%unitconv
        forall(i=1:limax,j=1:ljmax) &
          hourly(i,j)=(z_bnd(i,i,ik)-z_bnd(i,i,ik+1))*unit_conv

      case("COLUMN")    ! Column output in ug/m2, ugX/m2, molec/cm2
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        if(ih>0) hourly(:,:) = 0.0
        do iik=ik,KMAX_MID
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)                     &
                      + xn_adv(ispec,i,j,iik)             &
                      *                  unit_conv        & ! Units conv.
                      * roa(i,j,iik,1)                    & ! density.
                      * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
          endforall
        enddo

      case("COLUMNgroup")! GROUP Column output in ug/m2, ugX/m2, molec/cm2
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        if(ih>1) hourly(:,:) = 0.0
        do iik=ik,KMAX_MID
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)                     &
                      + dot_product(xn_adv(gspec,i,j,iik),&
                                    gunit_conv(:))        & ! Units conv.
                      * roa(i,j,iik,1)                    & ! density.
                      * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
          endforall
        enddo
        if(DEBUG) &
          write(*,'(a10,i7,a10,i7)')"K-level", ik, trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("SHLmcm3")        ! No cfac for short-lived species
        itot = ispec
        name = species(itot)%name
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                      * hr_out(ih)%unitconv  ! Units conv.
        endforall

      case("T2_C")        ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = t2_nwp(i,j,1) - 273.15

      case("rh2m")        ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = rh2m(i,j,1)*100

      case("ws_10m")      ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = ws_10m(i,j,1)
   
      case("heatsum")
       !hr_out(ih)%unit='degree_days'
        ik = hr_out(ih)%spec
        if(allocated(heatsum))then
          call CheckStop(ik,[1,size(heatsum,DIM=3)],"Hourly_out: '"//&
            trim(hr_out(ih)%type)//"' out of bounds!")
          forall(i=1:limax,j=1:ljmax) hourly(i,j)=heatsum(i,j,ik)
        else
          hourly(:,:) = 0.0
        endif

      case("pollen_left")
       !hr_out(ih)%unit='grains/m3'
        ik = hr_out(ih)%spec
        if(allocated(pollen_left))then
          call CheckStop(ik,[1,size(pollen_left,DIM=3)],"Hourly_out: '"//&
            trim(hr_out(ih)%type)//"' out of bounds!")
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = pollen_left(i,j,ik)
        else
          hourly(:,:) = 0.0
        endif

      case("pollen_emiss")
       !hr_out(ih)%unit='grains/m2/h'
        ik = hr_out(ih)%spec
        if(allocated(AreaPOLL))then
          call CheckStop(ik,[1,size(AreaPOLL,DIM=3)],"Hourly_out: '"//&
            trim(hr_out(ih)%type)//"' out of bounds!")
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = AreaPOLL(i,j,ik)
        else
          hourly(:,:) = 0.0
        endif

      case("theta")       ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = th(i,j,KMAX_MID,1)

      case("PRECIP")      ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = surface_precip(i,j)

      case("Idirect")     ! Direct radiation (W/m2); Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = Idirect(i,j)

      case("Idiffus")     ! Diffuse radiation (W/m2); Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = Idiffuse(i,j)

      case("D2D_inst")
       ! Here ispec is the index in the f_2d arrays
        call CheckStop(ispec<1.or.ispec>num_deriv2d,&
          "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%name)//", wrong D2D id!")
        if(hr_out(ih)%unit=="") hr_out(ih)%unit = f_2d(ispec)%unit
        unit_conv = hr_out(ih)%unitconv*f_2d(ispec)%scale

        if(debug_flag) then
          i=debug_li
          j=debug_lj
          write(*,"(2a,2i4,a,3g12.3)") "OUTHOUR "//trim(hr_out_type),&
            trim(hr_out(ih)%name), ih, ispec,trim(f_2d(ispec)%name),&
            d_2d(ispec,i,j,IOU_YEAR), d_2d(ispec,i,j,IOU_YEAR_LASTHH),&
            unit_conv
        endif
        forall(i=1:limax,j=1:ljmax)
          hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * unit_conv
        endforall
        if(debug_flag) &
          write(*,'(a,2i3,2es12.3)')"HHH DEBUG D2D", ispec, ih, &
            hr_out(ih)%unitconv, hourly(debug_li,debug_lj)

      case("D2D_accum")
       ! Here ispec is the index in the f_2d arrays
        call CheckStop(ispec<1.or.ispec>num_deriv2d,&
          "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%name)//", wrong D2D id!")
        if(hr_out(ih)%unit=="") hr_out(ih)%unit = f_2d(ispec)%unit
        unit_conv = hr_out(ih)%unitconv*f_2d(ispec)%scale
  
        if(debug_flag) then
          i=debug_li
          j=debug_lj
          write(*,"(2a,2i4,a,3g12.3)") "OUTHOUR "//trim(hr_out_type),&
            trim(hr_out(ih)%name), ih, ispec,trim(f_2d(ispec)%name),&
            d_2d(ispec,i,j,IOU_YEAR), d_2d(ispec,i,j,IOU_YEAR_LASTHH),&
            unit_conv
        endif
        forall(i=1:limax,j=1:ljmax)
          hourly(i,j) = (d_2d(ispec,i,j,IOU_YEAR)&
                        -d_2d(ispec,i,j,IOU_YEAR_LASTHH)) * unit_conv
          d_2d(ispec,i,j,IOU_YEAR_LASTHH)=d_2d(ispec,i,j,IOU_YEAR)
        endforall
        if(debug_flag) &
          write(*,'(a,2i3,2es12.3)')"HHH DEBUG D2D", ispec, ih, &
            hr_out(ih)%unitconv, hourly(debug_li,debug_lj)

      case("D2D_mean")
       ! Here ispec is the index in the f_2d arrays
        call CheckStop(ispec<1.or.ispec>num_deriv2d,&
          "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%name)//", wrong D2D id!")
        if(hr_out(ih)%unit=="") hr_out(ih)%unit = f_2d(ispec)%unit
        unit_conv = hr_out(ih)%unitconv*f_2d(ispec)%scale &
                   /MAX(nav_2d(ispec,IOU_YEAR)-navg(ih),1)
        navg(ih)=nav_2d(ispec,IOU_YEAR)
        
        if(debug_flag) then
          if(f_2d(ispec)%avg)then           ! averaged variables        
            write(*,"(a,1x,$)") "OUTHOUR D2D_mean avg"
          else                              ! accumulated variables --> mean
            write(*,"(a,1x,$)") "OUTHOUR D2D_mean acc"
          endif
          i=debug_li
          j=debug_lj
          write(*,"(a,2i4,a,3g12.3)"),&
            trim(hr_out(ih)%name), ih, ispec,trim(f_2d(ispec)%name),&
            d_2d(ispec,i,j,IOU_YEAR), d_2d(ispec,i,j,IOU_YEAR_LASTHH),&
            unit_conv
        endif
        forall(i=1:limax,j=1:ljmax)
          hourly(i,j) = (d_2d(ispec,i,j,IOU_YEAR)&
                        -d_2d(ispec,i,j,IOU_YEAR_LASTHH)) * unit_conv
          d_2d(ispec,i,j,IOU_YEAR_LASTHH)=d_2d(ispec,i,j,IOU_YEAR)
        endforall
        if(debug_flag) &
          write(*,'(a,2i3,2es12.3)')"HHH DEBUG D2D", ispec, ih, &
            hr_out(ih)%unitconv, hourly(debug_li,debug_lj)

      case DEFAULT
        call CheckStop( "ERROR-DEF! Hourly_out: '"//trim(hr_out(ih)%type)//&
                        "' hourly type not found!")

      endselect OPTIONS

      if(debug_flag) &
        write(*,"(a,3i4,2g12.3)")"DEBUG-HOURLY-OUT:"//trim(hr_out(ih)%name),&
           me,ih,ispec, hourly(debug_li,debug_lj), hr_out(ih)%unitconv

      !/ Get maximum value of hourly array
      hourly(limax+1:LIMAX,:) = 0.0
      hourly(1:limax,ljmax+1:LJMAX) = 0.0
      arrmax = maxval(hourly)

      if((hr_out(ih)%max>0.0).and.(arrmax>hr_out(ih)%max)) then
        write(*,*) "Hourly value too big!: ", ih, trim(hr_out(ih)%type), arrmax
        write(*,*) "Species : ", trim(name)," : ",  " ispec ", ispec
        write(*,*) "max allowed is : ",  hr_out(ih)%max
        write(*,*) "unitconv was   : ", hr_out(ih)%unitconv
        write(*,*) " me, limax, ljmax, LIMAX,LJMAX : ",  me, &
                          limax, ljmax ,LIMAX,LJMAX
        maxpos = maxloc(hourly)
        write(*,*) "Location is i=", maxpos(1), " j=", maxpos(2)
        write(*,*) "EMEP coords ix=", i_fdom(maxpos(1)), " iy=", j_fdom(maxpos(2))
        write(*,*) "hourly is ", hourly(maxpos(1),maxpos(2))
        if(hr_out(ih)%type(1:3)=="ADV") then
          write(*,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
          write(*,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
        endif
        errmsg="Error, Output_hourly/hourly_out: too big!"
      endif

!NetCDF hourly output
      def1%name=hr_out(ih)%name
      def1%unit=hr_out(ih)%unit
      def1%class=hr_out(ih)%type
      ist = max(IRUNBEG,hr_out(ih)%ix1)
      jst = max(JRUNBEG,hr_out(ih)%iy1)
      ien = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
      jen = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
      nk  = min(KMAX_MID,hr_out_nk)
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.

      select case(nk)
      case(1)       ! write as 2D
        call Out_netCDF(IOU_HOUR,def1,2,1,hourly,scale,CDFtype,ist,jst,ien,jen,ncFileID_given=ncFileID)
      case(2:)      ! write as 3D
        klevel=ik
        if(nk<KMAX_MID)  klevel=KMAX_MID-ik+1 !count from ground and up
        if(SELECT_LEVELS_HOURLY) klevel=k     !order is defined in LEVELS_HOURLY
        call Out_netCDF(IOU_HOUR,def1,3,1,hourly,scale,CDFtype,ist,jst,ien,jen,klevel,ncFileID_given=ncFileID)
      !case default   ! no output
      endselect
    enddo KVLOOP    
  enddo HLOOP
  ! CF convention: surface pressure to define vertical coordinates.
  if(NHOURLY_OUT>0.and.NLEVELS_HOURLY>0.and.all(hr_out(:)%name/="PS"))then
    def1%name='PS'
    def1%unit='hPa'
    def1%class='Surface pressure'
    CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
    scale=1.
    call Out_netCDF(IOU_HOUR,def1,2,1,ps(:,:,1)*0.01,scale,CDFtype,ist,jst,ien,jen,ncFileID_given=ncFileID)     
  endif

!Not closing seems to give a segmentation fault when opening the daily file
!Probably just a bug in the netcdf4/hdf5 library.
  call CloseNetCDF
  call CheckStop(errmsg)

! Write text file wich contents
  if(.not.(FORECAST.and.MasterProc))return
  i=index(filename,'.nc')-1;if(i<1)i=len_trim(filename)
  open(IO_HOURLY,file=filename(1:i)//'.msg',position='append')
  write(IO_HOURLY,*)date2string('FFFF: YYYY-MM-DD hh',current_date)
  close(IO_HOURLY)
endsubroutine hourly_out
