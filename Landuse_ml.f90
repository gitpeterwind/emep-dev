!> <Landuse_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
!! ************************************************************************! 

module Landuse_ml

use CheckStop_ml,   only: CheckStop,StopAll
use DO3SE_ml,       only: fPhenology, Init_DO3SE
use emep_Config_mod, only : LandCoverInputs 
use GridAllocate_ml,only: GridAllocate
use GridValues_ml,  only:  glat , glon   & ! latitude,
                          , i_fdom, j_fdom   & ! coordinates
                          , i_local, j_local &
                          , debug_proc, debug_li, debug_lj
use Io_ml,          only: open_file, ios, Read_Headers, Read2DN, IO_TMP &
                         ,IO_DO3SE
use KeyValueTypes,    only: KeyVal,KeyValue, LENKEYVAL
use LandDefs_ml,    only: Init_LandDefs, LandType, LandDefs, &
                          STUBBLE, Growing_Season,&
                          NLANDUSE_EMEP
use LandPFT_ml,       only: MapPFT_LAI, pft_lai
use ModelConstants_ml,only: DEBUG, NLANDUSEMAX, &
                            ! PALEO_TEST, & 
                            SEA_LIMIT, & 
                            USES, emep_debug, &
                            FLUX_VEGS,  nFluxVegs, & 
                            VEG_2dGS, VEG_2dGS_Params, & 
                            NPROC, IIFULLDOM, JJFULLDOM, &
                            DomainName, MasterProc
use MPI_Groups_ml, only : MPI_INTEGER,MPI_COMM_CALC, IERROR
use Par_ml,         only: LIMAX, LJMAX, &
                          limax, ljmax, me
!use Paleo_ml,       only: SetPaleo
use SmallUtils_ml,  only: wordsplit, find_index, NOT_FOUND, WriteArray, trims
use TimeDate_ml,    only: effectivdaynumber, nydays, current_date

use netcdf
use NetCDF_ml, only  : ReadField_CDF,check,printcdf

implicit none
private


!/- subroutines:

  public :: InitLanduse
  public :: ReadLanduse
  public :: SetLanduse
  private :: Polygon         ! Used for LAI
  private :: MedLAI          ! Used for LAI, Medit.

 integer, public, parameter :: NLUMAX = 30 ! max no. landuse per grid
 integer, private, save :: NLand_codes = 0 ! no. landuse in input files

! The LC: entries in the netcdf landuse, or the 
! the headers read from Inputs.Landuse define the "master-list" of
! codes for landuse. Each code must be present in the subsequent
! data files for phenology and DO3SE.

 character(len=20), dimension(NLANDUSEMAX), &
           public, save :: Land_codes = " " ! As used

 !=============================================
 type, public :: LandCov
   integer                   :: ncodes     ! Number of codes in grid
   integer,dimension(NLUMAX) :: &
          codes     &! landcover codes
         ,SGS       &! Start of growing season (days)
         ,EGS       &! End of growing season (days)
         ,Anth      &! Anthesis. For wheat
         ,Astart    &! Start photosynthetic activity, for DO3SE
         ,Aend       ! 
   real,   dimension(NLUMAX) :: &
          fraction  &! (coverage)
         ,LAI       &! Leaf-area-index (m2/m2)
         ,SAI       &! Surface-area-index (m2/m2) (leaves+bark, etc.)
         ,hveg      &! Max. height of veg.
         ,fphen     &! Potential (age) factor for Jarvis-calc
         ,Eiso      &! Emission potential, isoprene
         ,Emt       &! Emission potential, monoterpenes
         ,SumVPD    &! For critical VPD calcs, reset each day
         ,old_gsun   ! also for flux
 end type LandCov
 !=============================================
 type(LandCov), public, save, allocatable,dimension(:,:) :: LandCover
 !=============================================


 logical, public,save, allocatable,dimension(:,:) :: likely_coastal 
 logical, public,save, allocatable,dimension(:,:) :: mainly_sea 

 real,public,save, allocatable,dimension(:,:) :: water_fraction, ice_landcover 
 logical,public,save :: water_frac_set = .false.

 character(len=80), private :: errmsg

! integer ::ierror,mpi_comm_calc, mpi_integer
contains

 !==========================================================================
  subroutine InitLanduse(daynumber)
    integer, intent(in) :: daynumber
    logical :: filefound
    integer ::i,j,ilu,lu, ipar
    logical :: debug_flag = .false.
    ! Some config options, for veg with 2-D growing seasons
    ! needs map2d.. as real for ReadCDF
    integer :: n2dGS, n2dGSpars, i2dGS
    real,dimension(:,:,:),allocatable :: map2dGrowingSeasons
    character(len=len(VEG_2dGS_Params(1))) :: fname
    character(len=20) :: varname
    character(len=*), parameter :: dtxt='InitLanduse:'
    !=====================================

    !ALLOCATE ARRAYS
    allocate(LandCover(LIMAX,LJMAX))
    allocate(likely_coastal(LIMAX,LJMAX) )
    allocate(mainly_sea(LIMAX,LJMAX) )
    allocate(water_fraction(LIMAX,LJMAX), ice_landcover(LIMAX,LJMAX))

    ! First, check the number of "extra" (fake) vegetation 
    nFluxVegs = 0
    do ilu = 1, size( FLUX_VEGS )
        if(len_trim(FLUX_VEGS(ilu))>0) nFluxVegs=nFluxVegs+1
    end do

    if(MasterProc) write(*,*) dtxt//" nFluxVegs= ",nFluxVegs

    !ReadLandUse_CDF to be used as default when glc2000 data is improved?

    filefound=.false.
    call ReadLandUse(filefound) !=> Land_codes, Percentage cover per grid

    !ReadLandUse_CDF use Max Posch 5km landuse over emep area and glc200 where this dat is not defined.
    if(.not.filefound) then
        if(MasterProc) write(*,*) dtxt//" Into CDF "
        call ReadLandUse_CDF(filefound) !=> Land_codes, Percentage cover per grid
    end if


    ! Quick safety check
    ! we check that the length of the land-codes isn't equal to our declared
    ! length. That is usualyl a sign that the codes are too long and we may
    ! get into truncation worries.
    call CheckStop(maxval( len_trim(Land_codes(:))) >= len(Land_codes(1)),&
          "Land_Codes: increase size of character array" )

    call CheckStop(.not.filefound,"InitLanduse failed!")

    if(MasterProc) then
        print *,  dtxt//" Into Init_LandDefs ", NLand_codes
        print *,  dtxt//" Codes: ", Land_codes
        print *,  dtxt//" LandCoverInputs: ", LandCoverInputs
    end if
    !call Init_LandDefs(NLand_codes, Land_codes)   ! => LandType, LandDefs
    call Init_LandDefs(LandCoverInputs%LandDefs,NLand_codes, Land_codes)   ! => LandType, LandDefs

    !------ 2D maps of growing season, if set in config ----------------------- 

     n2dGS    = count( VEG_2dGS(:) /= "" )

    if ( n2dGS > 0 )  then

       n2dGSpars = count( VEG_2dGS_Params(:) /= "" ) - 1 ! First is fname

       allocate(map2dGrowingSeasons(LIMAX,LJMAX,n2dGS*n2dGSpars))

       fname = VEG_2dGS_Params(1)
       if(MasterProc) write(*,*) "CRU GS MAP ", trim(fname), n2dGS, n2dGSpars

       i2dGS = 0 
       do ilu = 1, n2dGS
         do ipar = 2, n2dGSpars+1

           varname = trims (&
               VEG_2dGS(ilu) // "_" // VEG_2dGS_Params(ipar) // ".txt" )

           i2dGS = i2dGS + 1 

           call ReadField_CDF(fname,varname, map2dGrowingSeasons(1,1,i2dGS),1, &
              interpol='zero_order',needed=.true.,debug_flag=.false.)  ! UnDef??

           if( debug_proc ) write(*,"(a20,5i5)") '2dGS '//trim(varname), ilu,&
              debug_li,debug_lj, i2dGS, &
               nint( map2dGrowingSeasons(debug_li,debug_lj,i2dGS) )

         end do ! ipar
       end do ! ilu
    end if
    !------  End 2D maps of growing season ------------------------------------- 


    ! effectiv daynumber to shift 6 month when in southern hemisphere
    effectivdaynumber=daynumber


    !/ -- Calculate growing seasons where needed and water_fraction
    !          (for Rn emissions)

    water_fraction(:,:) = 0.0
    ice_landcover(:,:)  = 0.0  !for Pb210 
    likely_coastal(:,:) = .false.
    mainly_sea(:,:)     = .false.

    do i = 1, limax
       do j = 1, ljmax

          debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj )
          do ilu= 1, LandCover(i,j)%ncodes
             lu      = LandCover(i,j)%codes(ilu)
             call CheckStop( lu < 0 .or. lu > NLANDUSEMAX , &
                  "SetLandUse out of range" )

             if ( LandDefs(lu)%SGS50 > 0 ) then ! need to set growing seasons 

                call Growing_season( lu,abs(glat(i,j)),&  
                     LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu) )
             else
                LandCover(i,j)%SGS(ilu) =  LandDefs(lu)%SGS50
                LandCover(i,j)%EGS(ilu) =  LandDefs(lu)%EGS50
             end if
            ! Hard-coded for now
             LandCover(i,j)%ANTH(ilu) =  -1
             if( LandDefs(lu)%name == "WinterWheat" ) then
                LandCover(i,j)%SGS(ilu) =  nint( map2dGrowingSeasons(i,j,1) )
                LandCover(i,j)%ANTH(ilu) =  nint( map2dGrowingSeasons(i,j,2) )
                LandCover(i,j)%EGS(ilu) =  nint( map2dGrowingSeasons(i,j,3) )
             else if( LandDefs(lu)%name == "SpringWheat" ) then
                LandCover(i,j)%SGS(ilu) =  nint( map2dGrowingSeasons(i,j,4) )
                LandCover(i,j)%ANTH(ilu) =  nint( map2dGrowingSeasons(i,j,5) )
                LandCover(i,j)%EGS(ilu) =  nint( map2dGrowingSeasons(i,j,6) )
!FUSK!
                LandCover(i,j)%SGS(ilu) =  LandCover(i,j)%EGS(ilu) - 90
             end if

             if ( DEBUG%LANDUSE>0 .and. debug_flag ) &
                  write(*,"(a,i3,a20,3i4)")"LANDUSE: LU_SETGS", &
                  lu, LandDefs(lu)%name,&
                  LandCover(i,j)%SGS(ilu),LandCover(i,j)%ANTH(ilu), &
                  LandCover(i,j)%EGS(ilu)


             !/ for landuse classes with bulk-resistances, we only
             !  need to specify height once. Dummy values are assigned
             !  to LAI and gpot:

             if ( LandType(lu)%is_bulk ) then
                LandCover(i,j)%hveg(ilu) =  LandDefs(lu)%hveg_max
                LandCover(i,j)%LAI(ilu)  =  0.0          
                LandCover(i,j)%fphen(ilu) =  0.0          
             end if

             if ( LandType(lu)%is_water ) water_fraction(i,j) = &
                  LandCover(i,j)%fraction(ilu)
             if ( LandType(lu)%is_ice   ) ice_landcover(i,j) = &
                  LandCover(i,j)%fraction(ilu)


          end do ! ilu

         ! Typically, we define as mainly sea when > 50% water, and
         ! likely_coastal when > 20%. SEA_LIMIT stores these numbers
         ! (We don't want to trust some squares with a mixture of sea
         !  and land for micromet purposes, e.g. T2 can be very wrong
         !  We mark these as likely coastal.)
         ! Unfortunately, we cannot yet determine if true sea or water

          if(water_fraction(i,j)>SEA_LIMIT(2) ) mainly_sea(i,j) = .true.
          if(water_fraction(i,j)>SEA_LIMIT(1).and. &
             water_fraction(i,j) < 0.999 ) likely_coastal(i,j) = .true.

          if ( DEBUG%LANDUSE>0 .and. debug_flag )  then
             write(*,"(a,2i4,f7.3,2L2)") "SEACOAST ", i_fdom(i), j_fdom(j), &
                water_fraction(i,j), mainly_sea(i,j), likely_coastal(i,j)
          end if
       end do ! j
    end do ! i

    water_frac_set = .true.  ! just to inform other routines
    
  end subroutine InitLanduse
 !==========================================================================
  subroutine ReadLanduse(filefound)

   logical :: filefound
   integer :: i,j,lu, index_lu, maxlufound
   character(len=20), dimension(NLANDUSEMAX+10) :: Headers
   type(KeyVal), dimension(10)      :: KeyValues ! Info on units, coords, etc.
   character(len=50) :: fname
   integer :: NHeaders, NKeys, Nlines
   logical :: debug_flag
   real :: sumfrac
   
  ! Specify the assumed coords and units - Read2DN will check that the data
  ! conform to these.
    type(keyval), dimension(2) :: CheckValues = &
        (/ keyval("Units","PercentGrid"), &
           keyval("Coords","ModelCoords") /)

 ! temporary arrays used.  Will re-write one day....
   real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_tmp ! tmp, with all data
   real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_in ! tmp, with all data
   real, dimension(LIMAX,LJMAX,NLUMAX):: landuse_data ! tmp, with all data
   integer, dimension(LIMAX,LJMAX):: landuse_ncodes ! tmp, with all data
   integer, dimension(LIMAX,LJMAX,NLUMAX):: landuse_codes ! tmp, with all data

   if ( DEBUG%LANDUSE>0 .and. MasterProc ) &
        write(*,*) "LANDUSE: Starting ReadLandUse "

   maxlufound = 0   
   Nlines = 0

   landuse_ncodes(:,:)   = 0     !***  initialise  ***
   landuse_codes(:,:,:)  = 0     !***  initialise  ***
   landuse_data  (:,:,:) = 0.0   !***  initialise  ***

!------------------------------------------------------------------------------

      ! Read Header info - this will define landuse classes for model

      fname = "Inputs.Landuse"
      if ( MasterProc ) then
         call open_file(IO_TMP,"r",fname,needed=.false.)
      end if
      call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_CALC,IERROR)
      if(ios==0)then
         if ( DEBUG%LANDUSE>0 .and. MasterProc ) write(*,*)'found '//trim(fname) 
         filefound=.true.

         call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,&
                     Headers,Keyvalues,CheckValues)
         
         call CheckStop( errmsg , "Read Headers" // fname )
         
         ! The first two columns are assumed for now to be ix,iy, hence:

         NHeaders = NHeaders -2
         call CheckStop( NHeaders /= NLANDUSE_EMEP, &
              "Inputs.Landuse not consistent with NLANDUSE_EMEP")

          NLand_codes=NHeaders        
        
         ! *** HERE we set the Landuse_codes ***
         do i = 1,  NLand_codes
            Land_codes(i) = trim ( Headers(i+2) )
         end do
         if(MasterProc)write(*,*)NLand_codes,' landuse categories defined from Inputs.Landuse:'
         if(MasterProc)write(*,fmt="(20(A,1x))")(trim(Land_codes(i)),i=1,NLand_codes)

         ! Then data:
         
         call Read2DN("Inputs.Landuse",NLand_codes,landuse_in,&
                 HeadersRead=.true.)
         
         !-------------------------------------------------------------------
         
         if ( DEBUG%LANDUSE>0 .and. MasterProc ) then
            write(*,*) "LANDUSE: LAND_CODES ARE ", NHeaders
            call WriteArray(Land_codes,NLand_codes,"Land_Codes")
         end if
         
      else
         filefound=.false.
         if(MasterProc)Write(*,*)'Inputs.Landuse not found'
         return
         call StopAll('Inputs.Landuse not found') 
      end if

!      call printCDF('LU', landuse_in(:,:,1),'??')

    do i = 1, limax
       do j = 1, ljmax
           debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
           do lu = 1, NLand_codes
              if ( landuse_in(i,j,lu) > 0.0 ) then

                 call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                         index_lu, maxlufound, landuse_codes, landuse_ncodes)
   
                     landuse_data(i,j,index_lu) = &
                       landuse_data(i,j,index_lu) + 0.01 * landuse_in(i,j,lu)
               end if
               if ( DEBUG%LANDUSE>0 .and. debug_flag )  &
                       write(*,"(a15,i3,f8.4,a10,i3,f8.4)") "DEBUG Landuse ",&
                          lu, landuse_in(i,j,lu), &
                           "index_lu ", index_lu, landuse_data(i,j,index_lu)
           end do ! lu
           LandCover(i,j)%ncodes  = landuse_ncodes(i,j)
           LandCover(i,j)%codes(:) = landuse_codes(i,j,:)
           LandCover(i,j)%fraction(:)  = landuse_data(i,j,:)

           sumfrac = sum( LandCover(i,j)%fraction(:) )

             if (  sumfrac < 0.99 .or. sumfrac > 1.01 ) then
               write(unit=errmsg,fmt="(a19,3i4,f12.4,8i4)") &
                 "Land SumFrac Error ", me,  &
                    i_fdom(i),j_fdom(j), sumfrac, limax,  ljmax, &
                       i_fdom(1), j_fdom(1), i_fdom(limax), j_fdom(ljmax)
               call CheckStop(errmsg)
             end if

      end do  !j
   end do  !i

   if (DEBUG%LANDUSE>1) write(6,*) "Landuse_ml: me, Nlines, maxlufound, ascii = ", &
                                  me, Nlines, maxlufound

  end subroutine  ReadLanduse
 
  subroutine ReadLanduse_CDF(filefound)
    !Read data in other grid and interpolate to present grid
    !
    !So far only basic version for use in TNO7. Under construction
    !
    implicit none
    logical :: filefound
    integer :: i,j,lu, ilu, index_lu, maxlufound, iam, iveg, ifile
    logical :: debug_flag
    real :: sumfrac

    character(len=40) :: varname
    character(len=200) :: fname1,fname2, fName
!XX    character(len=200), dimension(2) :: fnamesLC
    integer :: ncFileID, nDimensions,nVariables,nAttributes,timeDimID,varid
    integer :: nwords, err, xtype,ndims  ,status
    character(len=20) :: ewords(7), code ! LC:CF:EMEP
    character(len=*), parameter :: dtxt='RdLanduseCDF:'
    logical :: fexist=.false.!file exist flag
  
    real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_in ! tmp, with all data
    real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_glob  ! CLM crude
    real, dimension(LIMAX,LJMAX):: landuse_tot ! CLM 
    real, dimension(LIMAX,LJMAX):: landuse_tmp ! tmp, with all data
    !CLM cude:
    logical :: dbgclm
    integer :: mdbg, idbg, jdbg  ! debug i,j
    real    :: dbgsum
    logical :: EuroFileFound=.false., GlobFileFound=.false. !file exist flags

    real, dimension(LIMAX,LJMAX,NLUMAX):: landuse_data ! tmp, with all data
    integer, dimension(LIMAX,LJMAX):: landuse_ncodes ! tmp, with all data
    integer, dimension(LIMAX,LJMAX,NLUMAX):: landuse_codes ! tmp, with all data
    logical, save :: dbg0=.false.

    if(  DEBUG%LANDUSE>0 .and. MasterProc )  then
       write(*,*) dtxt//" Starting"
       dbg0 = .true.
    end if
          mdbg = 10; idbg = 1; jdbg = 1
          mdbg =  9; idbg = 5; jdbg = 11
          dbgclm= ( me==mdbg ) !!! .and. i_fdom(i)==90 .and. j_fdom(j)==106) 


    if (MasterProc ) write(*,*) dtxt//"LANDUSE_CDF:"
    !    filefound=.false.
    !    return

    maxlufound = 0   

    landuse_ncodes(:,:)   = 0     !***  initialise  ***
    landuse_codes(:,:,:)  = 0     !***  initialise  ***
    landuse_data  (:,:,:) = 0.0   !***  initialise  ***
    landuse_in  = 0.0              !***  initialise  ***
    landuse_glob  = 0.0              !***  initialise  ***

    !Landusefile where landcodes are not predefined, but read from the file.
    !DS fName1='Landuse_PS_5km_LC.nc'
    !DS fName2='LanduseGLC.nc'
    !1)check that file exists
    !note that every processor open and read the same file
    ilu=0 ! DS here
    !XX fnamesLC = [ 'Landuse_PS_5km_LC.nc', 'LanduseGLC.nc' ]

    if(MasterProc) print *,'LC NAMES ', &
      trim(LandCoverInputs%MapFile(1)), &
      trim(LandCoverInputs%MapFile(2)) 

    !FILELOOP: do ifile = 1, size(fnamesLC)
    FILELOOP: do ifile = 1, size(LandCoverInputs%MapFile(:))
      fName = LandCoverInputs%MapFile(ifile)   !XX fnamesLC(ifile)
      status=nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID)
      inquire(file=trim(fName),exist=fexist)
      if ( dbg0 .and. fexist)write(*,*) dtxt//"LANDUSE: found "//trim(fName)

      if(status /= nf90_noerr) then
        !XX print *, dtxt//"XXXXXX", ifile, trim(fnamesLC(ifile)), trim(fName)
         if ( dbg0 )write(*,*) dtxt//"LANDUSE: NOT found "//trim(fName)
         cycle !CLM
      end if
      if ( ifile==1 ) EuroFileFound = .true. ! CRUDE; DANGEROUS; ARGH!
      if ( ifile==2 ) GlobFileFound = .true. ! CRUDE; DANGEROUS; ARGH!
      if ( dbg0 )write(*,*) dtxt//"LANDUSE: found "//trim(fName)

     !get list of variables
      call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,timeDimID))
      ! All the inquire functions are inexpensive to use and require no I/O, since the information
      ! they provide is stored in memory when a netCDF dataset is first opened.   

      !loop over all variables in file
      VARIDLOOP1: do varid=1,nVariables
          if ( DEBUG%LANDUSE>0 )  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

          call check(nf90_Inquire_Variable(ncFileID,varid,varname,xtype,ndims))
          if ( dbg0 )write(*,*) dtxt//"checking "//trim(fName)//':'// &
                trim(varname), index( varname, "LC:") 

          ! landcover terms look like, e.g. LC:CF:EMEP
          if( index( varname, "LC:") < 1 ) then
     !       print *, dtxt//"Skips ", ifile, trim(varname)
            cycle ! ONLY LC: (LandCode) wanted
          end if
          call wordsplit(varname,3,ewords,nwords,err,separator=":")
          if( ewords(3) /= "EMEP" .and.  ewords(3) /= "CLM" ) cycle ! ONLY EMEP coded for now
     !       print *, dtxt//"Does ", me, ifile, trim(varname), trim(ewords(2))

          !=========================
          if( ewords(2) == "IAM_VEG" .and. nFluxVegs < 1 ) exit  ! No IAM veg to process
          !=========================

          !DS ilu=ilu+1
          !DS if ( dbg0 )&
          !DS      write(*,*) "defining new LC "//ewords(2)//"  ilu= " , ilu

!CHECK HERE to see if we already have this landcode:...
          lu = find_index( ewords(2), Land_codes(:) ) 
          call CheckStop( ilu>NLANDUSEMAX , &
               dtxt//"NLANDUSEMAX smaller than number of landuses defined in file "//trim(fname) )
          if(me==mdbg) print *, dtxt//'LUVARNAME', ifile, ilu, trim(varname)
          if(me==mdbg) print *, dtxt//'LULU', ifile, ilu, trim(ewords(1)), trim(ewords(2))
          if (  lu > 0 ) then
            if ( dbg0 ) write(*,*) dtxt//"Already have"//ewords(2), lu
          else
            ilu = ilu + 1
            lu  = ilu  
            if ( dbg0 ) write(*,*) dtxt//"Adding code"//ewords(2), ilu
            Land_codes(ilu) = ewords(2)    ! Landuse code found on file
          end if

          !TMP call ReadField_CDF(trim(fName),varname,& 
          !TMP      !CLM landuse_in(1,1,ilu),1,interpol='conservative', &
          !TMP      landuse_tmp,1,interpol='conservative', &
           !TMP     needed=.true.,debug_flag=.true.,UnDef=-9.9E19) 
               !TMP needed=.true.,debug_flag=.false.,UnDef=-9.9E19) 

          if ( ifile == 1 ) then
          call ReadField_CDF(trim(fName),varname,& 
               !CLM landuse_in(1,1,ilu),1,interpol='conservative', &
               landuse_tmp,1,interpol='conservative', &
               needed=.true.,debug_flag=.true.,UnDef=-9.9E19) 
               landuse_in(:,:,lu) = landuse_tmp
               landuse_tot(:,:) = landuse_tot(:,:) + landuse_tmp
             if ( dbgclm ) then
               print *, "F1 ", landuse_in(idbg,jdbg,lu), landuse_tot(idbg,jdbg), lu, trim(ewords(2))
               print *, "F1 latlon ", glat(idbg,jdbg), glon(idbg,jdbg)
             end if
          else
          call ReadField_CDF(trim(fName),varname,& 
               !CLM landuse_in(1,1,ilu),1,interpol='conservative', &
               landuse_tmp,1,interpol='conservative', &
               known_projection='lon lat', &
               needed=.true.,debug_flag=.true.,UnDef=-9.9E19) 
               landuse_glob(:,:,lu) = landuse_tmp
             if ( dbgclm ) then
               print *, "F2 ", landuse_in(idbg,jdbg,lu), &
                 landuse_tot(idbg,jdbg), landuse_glob(idbg,jdbg,lu), trim(ewords(2))
             end if
          end if

         ! ---- TMP. Will sort out GLC file  another day
         ! Some "IAM" veg species can be defined for calculations of ozone
         ! fluxes. These are assigned very small land-area, using the mask
         ! which the IAM_VEG species gives.  
         ! We divive the area  by the nFluxVegs to keep the total area small

          if ( Land_codes(ilu) == "IAM_VEG" ) then 

             varname = "IAM_DF" ! good enough
             iveg = ilu
             forall ( i=1:limax,j=1:ljmax)
                landuse_in(i,j,ilu) = landuse_in(i,j,ilu) / real(nFluxVegs)
             end forall

           IAM_VEG: do iam = 1, size( FLUX_VEGS )
             if ( len_trim( FLUX_VEGS(iam) ) < 1 ) then

                if(MasterProc) write(*,*)"Landuse SKIPS IAM ", iam
                cycle IAM_VEG
             end if
              
             ilu = iveg-1 + iam ! first iam overwrites IAM_VEG name
             if(MasterProc) write(*,*)"Landuse EXTRA IAM ",&
                iam, ilu, FLUX_VEGS(iam)
                    
             Land_codes(ilu) = FLUX_VEGS(iam)
             forall ( i=1:limax,j=1:ljmax)
                landuse_in(i,j,ilu) = landuse_in(i,j,iveg)
             end forall
           end do IAM_VEG
          end if ! IAM_VEG
         ! ------

         end do VARIDLOOP1
         call check(nf90_close(ncFileID))!fname1
       end do FILELOOP ! DSLC

       NLand_codes=ilu !DSCLM now here

!DSQUERY          if(fexist .and. any(landuse_in(1:limax,1:ljmax,ilu)<-0.1))then
!DSQUERY             !complete missing data with data from second file
!DSQUERY             !name in second file may be defined differently
!DSQUERY             varname=Land_codes(ilu) !name such as CF (without LC: etc.)


!DS             call ReadField_CDF(trim(fName2),varname,&
!DS                  landuse_tmp,1,interpol='conservative', &
!DS                  needed=.false.,debug_flag=.false.)
                  !CLM needed=.true.,debug_flag=.false.)

       if  ( EuroFileFound .and. GlobFileFound ) then ! we need to merge
             do j = 1, ljmax
                do i = 1, limax
                   !landuse_tmp can be numerically larger than 1.0 (1E-15 larger). That made negative deposition for PB210
                   !DS if(landuse_in(i,j,ilu)<-0.1)landuse_in(i,j,ilu)=min(1.0,landuse_tmp(i,j))
             if ( dbgclm .and. i==idbg.and.j==jdbg ) then
               print *, "F3 ", landuse_tot(idbg,jdbg)
             end if
                   if(landuse_tot(i,j)< 0.99999 ) then
                      landuse_in(i,j,:)= 0.0  ! Will overwrite all PS stuff
                      dbgsum = 0.0
                      do ilu = 1, NLand_codes
                        landuse_in(i,j,ilu) = min(1.0, landuse_glob(i,j,ilu) )
                         dbgsum = dbgsum + landuse_in(i,j,ilu)
             if ( dbgclm .and. i==idbg.and.j==jdbg ) then
               print "(a,i3,3es15.6,1x,a)", "F4 ", ilu, landuse_in(idbg,jdbg,ilu), &
                  landuse_tot(idbg,jdbg), dbgsum, trim(Land_Codes(ilu))
             end if
                      end do
                   end if
                end do  !j
             end do  !i
       end if 
!DSCLM          end if

          if(MasterProc) write(*,*)"LandDefs DONE ", ilu, Land_codes(ilu), maxval( landuse_in ), minval(landuse_in)
!DSCLM       NLand_codes=ilu
       if(MasterProc) then
            write( *,*) "Number of landuse codes ", NLand_codes
            write( *,*) "LAND_CODES: ", Land_codes(1:NLand_codes)
       end if

!DS    else
!DS       !the landusefile with softcoded lancodes has not been found. Use "old" method 
!DS       if ( dbg0 )write(*,*) "LANDUSE: LC: not found "//trim(fName1)
!DS       call CheckStop("Landuse: No landcover files")

 !DS   end if !switch hardcoded/fileread lu definitions

    do i = 1, limax
       do j = 1, ljmax
          do lu = 1, NLand_codes
             if ( landuse_in(i,j,lu) > 0.0 ) then

                call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                     index_lu, maxlufound, landuse_codes, landuse_ncodes)
                landuse_data(i,j,index_lu) = &
                     landuse_data(i,j,index_lu) + landuse_in(i,j,lu)!already in fraction unit
             end if
          end do ! lu
          LandCover(i,j)%ncodes  = landuse_ncodes(i,j)
          LandCover(i,j)%codes(:) = landuse_codes(i,j,:)
          LandCover(i,j)%fraction(:)  = landuse_data(i,j,:)
          sumfrac = sum( LandCover(i,j)%fraction(:) )


            if (  sumfrac < 0.99 .or. sumfrac > 1.01 ) then
               write(unit=errmsg,fmt="(a19,5i4,f12.4,8i4)") &
                 dtxt//" SumFrac Error ", me,i,j,  &
                    i_fdom(i),j_fdom(j), sumfrac, limax,  ljmax, &
                       i_fdom(1), j_fdom(1), i_fdom(limax), j_fdom(ljmax)
               write(*,*)dtxt//'lat/lon: ',trim(errmsg),glat(i,j), glon(i,j)
               if(abs(sumfrac-1.0)<0.2.and.abs(glat(i,j))>89.0)then
                  write(*,*)'WARNING: ',trim(errmsg),sumfrac,glat(i,j)
               else
                   write(*,*)'lat/lon: ',trim(errmsg),glat(i,j), glon(i,j)
                 call CheckStop(errmsg)
               end if
             end if

       end do  !j
    end do  !i


    filefound=.true.
   if (DEBUG%LANDUSE>0) write(6,*) "Landuse_ml: me,  maxlufound, cdf = ", &
                                  me, maxlufound

  end subroutine ReadLanduse_CDF

  !=========================================================================
  subroutine  SetLandUse(daynumber, month)
    integer, intent(in) :: daynumber, month
    integer :: i,j,ilu,lu ! indices
    integer, save :: old_month = -1
    integer, save :: old_daynumber = -1
    logical, save :: my_first_call = .true.
    logical, save :: init_needed=.true. ! since my_first_call had some confusions..
    logical :: debug_flag = .false., debug_sgs
    real :: hveg, lat_factor
    real :: xSAIadd
    integer :: pft
    logical, save :: debugProc = .false.
    character (len=*), parameter :: dtxt='SetLandUse:'
    character (len=60) :: dnam !mainly for debug

! Treatment of growing seasons in the southern hemisphere:
!   all the static definitions (SGS,EGS...) refer to northern hemisphere, 
!   but the actual simulation dates are shifted by 6 months in the southern
!   hemisphere by using uses effectivdaynumber and 
!   mod(current_date%month+5,12)+1 in southern hemis



   !======================================================================
    if ( my_first_call ) then

        if ( DEBUG%LANDUSE>0 .and. debug_proc ) debugProc = .true.

       !read in data from file
        my_first_call   = .false.
        
        call InitLanduse(daynumber)

       ! The DO3SE params are needed for the call to fPhenology
      
        call Init_DO3SE(IO_DO3SE, &
              LandCoverInputs%Do3seDefs,NLand_codes, Land_codes, errmsg)
      !trim(LandCoverInputs%Do3sDefs), &
      !  call Init_DO3SE(IO_DO3SE,"Inputs_DO3SE.csv",NLand_codes, Land_codes, errmsg)
        call CheckStop(errmsg, "Reading DO3SE ")

    end if ! my_first_call
   !======================================================================

    if ( daynumber == old_daynumber ) then
        my_first_call = .false. ! PW
        return
    end if
    old_daynumber = daynumber

    if(MasterProc) write(*,*) dtxt//" day, pfts? ", daynumber, USES%PFT_MAPS
    if(debugProc ) write(*,"(a,5i5,L2)") dtxt//" debug me i j pft? ", me, &
         debug_li, debug_lj, limax, ljmax, USES%PFT_MAPS


   !Landcover data can be set either from simplified LPJ
   !PFTs, or from the "older" DO3SE inputs file

     if ( USES%PFT_MAPS ) then !- Check for LPJ-derived data -
         if (MasterProc) print *, dtxt//"New PFTMAPS ", month, old_month
         if ( month /= old_month ) then 
           call MapPFT_LAI( month )
         end if
     end if
    
!PALEO LANDUSE
!    if( PALEO_TEST ) then
!     call SetPaleo(daynumber, month)
!    end if



     do i = 1, limax
       do j = 1, ljmax

          debug_sgs = ( DEBUG%LANDUSE > 1 .and. &
            current_date%hour == 0  .and.  &
            glat(i,j) <   1.0 .and. glat(i,j) > -1.0 .and.  &
            glon(i,j) > -72.0 .and. glon(i,j) < -70.0 )


          effectivdaynumber=daynumber
         ! effectiv daynumber to shift 6 months when in southern hemisphere
          if(glat(i,j)<0.0)effectivdaynumber=mod(daynumber+182,nydays)+1 

          debug_flag = ( debugProc .and. i == debug_li .and. j == debug_lj ) 

          if ( debug_flag ) then
                 write(*,"(a,i3,9i6)") dtxt//" debug DATE ", &
                  LandCover(i,j)%ncodes, daynumber, current_date
          end if

          do ilu= 1, LandCover(i,j)%ncodes
             lu      = LandCover(i,j)%codes(ilu)
             pft     = LandType(lu)%pft
             dnam    = dtxt//trim(LandDefs(lu)%name)

             if ( debug_flag ) print *, trim(dnam)//" lu pft", lu, pft,&
               LandType(lu)%is_bulk, LandType(lu)%is_forest

             if ( LandType(lu)%is_bulk ) then
                LandCover(i,j)%LAI(ilu) = 0.0
                LandCover(i,j)%SAI(ilu) = 0.0
                cycle    
             end if!else Growing veg present:

            if ( LandDefs(lu)%name == "MED_OAK" .or.  &
                  LandDefs(lu)%name == "MED_PINE"   ) then

                LandCover(i,j)%LAI(ilu) = MedLAI(effectivdaynumber, &
                   100, 166, & ! Hard-code from Mapping Manual
                     LandDefs(lu)%LAImin, LandDefs(lu)%LAImax )
                if ( debug_flag ) then
                   write(*,"(a,3i4,3f8.3)") "MED_TREE "//&
                     trim(LandDefs(lu)%name), effectivdaynumber,&
                     LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                     LandDefs(lu)%LAImin, LandDefs(lu)%LAImax, &
                     LandCover(i,j)%LAI(ilu)
                end if

             else
                LandCover(i,j)%LAI(ilu) = Polygon(effectivdaynumber, &
                    0.0, LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                      LandCover(i,j)%SGS(ilu), LandDefs(lu)%SLAIlen, &
                         LandCover(i,j)%EGS(ilu), LandDefs(lu)%ELAIlen)
             end if

             LandCover(i,j)%fphen(ilu) = fPhenology( lu &
                ,effectivdaynumber &
                ,LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)&
                ,debug_flag )

             if ( debug_flag ) then
             !if (debug_sgs  ) then
               write(*,"(a,3i4,5f8.3)")trim(dnam)//" CHECK_VEG ",&
                 effectivdaynumber, &
                 LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                 LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                 LandCover(i,j)%LAI(ilu), LandCover(i,j)%fphen(ilu)
               write(*,"(a,L3)")trim(dnam)//' CHECK_PFT', USES%PFT_MAPS
            end if

            if ( USES%PFT_MAPS ) then
                if ( DEBUG%PFT_MAPS > 0 .and. debug_flag ) then
                     if ( pft > 0 ) then
                       write(*,"(a,i4,i6,2f8.3)") trim(dnam)//" PFTS COMP? ", &
                          daynumber, pft, LandCover(i,j)%LAI(ilu), &
                           pft_lai(i,j, pft)*LandDefs(lu)%LAImax
                     else
                       write(*,"(2a,i4,i6,2f8.3)") trim(dnam)//" PFTS COMP? ", &
                          daynumber, pft,&
                           LandCover(i,j)%LAI(ilu), -1.0
                     end if
                 end if
                 if ( pft > 0 ) then !PFT OVERWRITE!
                    LandCover(i,j)%LAI(ilu)= pft_lai(i,j, pft)*LandDefs(lu)%LAImax
               if(debug_flag) write(*,"(a,2i4,5f8.3)")dtxt//' CHECK_LAI', lu, pft, &
                      pft_lai(i,j, pft),LandDefs(lu)%LAImax
                    LandCover(i,j)%fphen(ilu)= 1.0  ! Skip fphen if using PFT
                    LandCover(i,j)%SGS(ilu)=  -999  ! Marker, since not used
                    LandCover(i,j)%EGS(ilu)=  -999  ! Marker, since not used
                 end if
                if (debug_sgs  )write(*,*) "ESGS in PFTMAPS"
            end if


             hveg = LandDefs(lu)%hveg_max   ! defaults
             xSAIadd = 0.0

             if (  LandType(lu)%is_crop ) then

               ! Note that IAM crops have SLAIlen=0, so are immediately
               ! given LAI=3.5, SAI=5.

                if ( effectivdaynumber < LandCover(i,j)%SGS(ilu) .or. &
                     effectivdaynumber > LandCover(i,j)%EGS(ilu)  ) then
                   hveg = STUBBLE
                   xSAIadd = 0.0
                else if ( effectivdaynumber < &
                     (LandCover(i,j)%SGS(ilu) + LandDefs(lu)%SLAIlen) ) then
                   hveg=  LandDefs(lu)%hveg_max * &
                     LandCover(i,j)%LAI(ilu) / LandDefs(lu)%LAImax
                   xSAIadd = ( 5.0/3.5 - 1.0) * LandCover(i,j)%LAI(ilu)
                else if ( effectivdaynumber < LandCover(i,j)%EGS(ilu) ) then
                   hveg = LandDefs(lu)%hveg_max  ! not needed?
                   xSAIadd = 1.5  ! Sensescent
                end if
                LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + xSAIadd

             ! end if ! crops


           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

             else if( LandType(lu)%is_forest ) then
               if ( glat(i,j) >= 60.0 ) then
                       lat_factor  = max(0.3, ( 1.0 - 0.05* (glat(i,j)-60.0)) )
                       hveg  = hveg *  lat_factor
                       LandCover(i,j)%LAI(ilu) = LandCover(i,j)%LAI(ilu)  * lat_factor
               end if
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + 1.0
             else if( LandType(lu)%is_seminat ) then !A2017 SNL
               if ( glat(i,j) >= 60.0 ) then
                       lat_factor  = max(0.3, ( 1.0 - 0.05* (glat(i,j)-60.0)) )
                       hveg  = hveg *  lat_factor
                       LandCover(i,j)%LAI(ilu) = LandCover(i,j)%LAI(ilu)  * lat_factor
               end if
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + 0.5 !  A2017SNL
             else
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)   !defaults
             end if

             LandCover(i,j)%hveg(ilu) =  hveg


if( debug_sgs .or. debug_flag  ) then
   write(*, "(a20,i4,2f7.1,3i5,f8.2)") trim(dnam)//":ESGS:",&
     lu, glat(i,j), glon(i,j), &
     daynumber, effectivdaynumber, Landcover(i,j)%SGS(ilu), Landcover(i,j)%LAI(ilu)
!end if
!             if (debug_sgs  ) then
               write(*,"(a,3i4,5f8.3)")trim(dnam)//"CHECK_VEGB:",&
                 effectivdaynumber, &
                 LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                 LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                 LandCover(i,j)%LAI(ilu), LandCover(i,j)%fphen(ilu)
end if ! debug_sgs
         end do ! lu
       end do ! j
    end do ! i

    my_first_call   = .false.

! --- print out for debug cell
    if ( DEBUG%LANDUSE>0.and.debug_proc ) then
       i=debug_li
       j=debug_lj

       do ilu= 1, LandCover(i,j)%ncodes
          lu      = LandCover(i,j)%codes(ilu)
          pft     = LandType(lu)%pft
          if ( LandType(lu)%is_bulk ) cycle    !else Growing veg present:

            write(*,"(a,i3,a16,i4,f7.2,3f8.3,2i5)") "LANDUSE Phen ", lu,&
             trim(LandDefs(lu)%name), daynumber, LandCover(i,j)%hveg(ilu),&
              LandCover(i,j)%SAI(ilu), LandCover(i,j)%LAI(ilu), &
              LandCover(i,j)%fphen(ilu), &
             LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)
       end do

   end if
 end subroutine  SetLandUse
! =====================================================================

!=======================================================================
function Polygon(jdayin,Ymin,Ystart,Ymax,Sday,LenS,Eday,LenE) &
result (Poly)
!=======================================================================

!     Calculates the value of a parameter Y with a polygon
!     distribution - currently LAI and g_pot

!            _____________       <- Ymax
!           /             \
!          /               \
!         /                 \
!        /                   \
!       |                     |  <- Ystart
!       |                     |
!       |                     |
!  ----------------------------- <- Ymin
!       S  S1            E1   E
!

!   Inputs
    integer, intent(in) :: jdayin     !day of year
    real, intent(in) ::    Ymin       !minimum value of Y
    real, intent(in) ::    Ystart     !value Y at start of growing season
    real, intent(in) ::    Ymax       !maximum value of Y
    integer, intent(in) ::    Sday    !start day (e.g. of growing season)
    integer, intent(in) ::    LenS    !length of Start period (S..S1 above)
    integer, intent(in) ::    Eday    !end day (e.g. of growing season)
    integer, intent(in) ::    LenE    !length of end period (E..E1 above)

!  Output:
    real ::   Poly  ! value at day jday

! Local
    integer :: jday ! day of year, after any co-ordinate change
    integer ::    S ! start day
    integer ::    E ! end day
    
    jday = jdayin
    E = Eday
    S = Sday

  ! Here we removed a lot of code associated with the leaf-age
  ! version of g_pot. 
       
    if ( jday  <  S .or. jday >  E ) then
       Poly = Ymin
       return
    end if


    if (jday <=  S+LenS  .and. LenS > 0 ) then

        Poly = (Ymax-Ystart) * (jday-S)/LenS  + Ystart 

    else if ( jday >=  E-LenE .and. LenE > 0.0 ) then   !d1.1 test for LenE

        Poly = (Ymax-Ystart) * (E-jday)/LenE + Ystart

    else

        Poly =Ymax

    end if
    

 end function Polygon

 !=======================================================================
function MedLAI(jday,LAIs,LAIe,LAImin,LAImax) result (LAI)
!=======================================================================

!     Calculates the value of LAI from the Mapping manual
!     functions for Mediteranean forests

!   Inputs
    integer, intent(in) :: jday     !day of year, after any co-ordinate change
    integer, intent(in) ::  LAIs,LAIe
    real, intent(in) ::     LAImin,LAImax

!  Output:
    real ::   LAI  ! value at day jday


    if( jday <= LAIs ) then
      LAI = 0.35*((real(LAIs-jday))/LAIs) + LAImin
    else if ( jday < (366.0-LAIe) ) then
      LAI = (LAImax-LAImin)*(real(jday-LAIs)/LAIs) + LAImin
    else
      LAI = (LAImax-( LAImin+0.35))*(real(366-jday)/LAIe) + LAImin+0.35
    end if


 end function MedLAI

 !=======================================================================

end module Landuse_ml
