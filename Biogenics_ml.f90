module Biogenics_ml
  !/-- reads in forest data, defines natural emissions  arrays
  !
  ! Usage: o3mod.f calls init_forests
  !
  ! The rest of the biogenic stuff is calculated in Setup_1d
  ! To Do:
  !   move fac from setup to here - apply to emforest. This keeps
  !   canopy_ecf O(1).
  !---------------------------------------------------------------------------
  use Par_ml   , only : MAXLIMAX,MAXLJMAX
  use My_Emis_ml       , only : NFORESTVOC
  implicit none
  private

  !/-- subroutines
  public ::  Forests_Init      ! was rforest
!hf u2:
  logical, private, parameter:: DEBUG = .false.
  integer, public, parameter :: NBIO = NFORESTVOC
  integer, public, save      :: BIO_ISOP, BIO_TERP

  real, public, save, dimension(NFORESTVOC,MAXLIMAX,MAXLJMAX) :: &
      emforest    & !  Gridded standard (30deg. C, full light) emissions
     ,emnat         !  Gridded std. emissions after scaling with density, etc.

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers
  !    - light correction v. crude just now, = 1 or 0. Fix later..

  real, public, save, dimension(NFORESTVOC,40) :: canopy_ecf  ! Canopy env. factors
                                                        
  !/-- DMS factors

  integer, public, parameter :: IQ_DMS = 35  ! code for DMS emissions

  !u4 real, public,    save :: dms_fact
  logical, public, save :: first_dms_read

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    subroutine Forests_Init()

  use Par_ml   , only : GIMAX,GJMAX &
                ,ISMBEG,JSMBEG &
                ,NPROC,me,MSG_READ1 &
                ,li0,li1,lj0,lj1
  use My_Emis_ml       , only : FORESTVOC
!hf u2  use My_Runmode_ml    , only :  DEBUG, stop_test 
  use GridValues_ml    , only : iclass,xm2
  use Io_ml            , only : IO_FORES, open_file, ios
!    Read forest data for natural VOC emissions
!-----------------------------------------------------------------------------
!ds 001011 - slightly modified to use free-format, some F90, and
!     to obtain terpene emissions. Many old lines deleted
!     cemfa for isoprene renamed to ecf_isop. ecf_terp added for terpene
!-----------------------------------------------------------------------------


    integer i, j, d, info, ii, jj, it, cat
    integer itmp, jtmp, i50, j50
    real tmp(6)
    real sumland

    real, dimension(6,GIMAX,GJMAX)       ::  gforest  ! Forest on global domain
    real, dimension(6,MAXLIMAX,MAXLJMAX) ::  forest   ! Forest on local domain
    real :: forsum(6), isopsum, terpsum
    real agts, agr, agtm, agct1, agct2, agcl &
               , agct, oak, decid, conif, spruc, crops, sitka &
               ,itk, fac
    integer :: inp_error

    do i = 1, NFORESTVOC
      if ( FORESTVOC(i) == "isoprene" ) then
            BIO_ISOP = i
      else if ( FORESTVOC(i) == "terpene " ) then
            BIO_TERP = i
      else
          call gc_abort(me,NPROC,"BIO_ERROR")
!hf u2    call stop_test(.true.,me,NPROC,99,"BIO_ERROR")
      end if
    end do

    inp_error = 0
    if (me ==  0) then

        gforest(:,:,:) = 0.0

      !ds.......................start of new 1.................................
      !     Read in forest Land-use data in percent

        call open_file(IO_FORES,"r","forest.pcnt",needed=.true.)
        if (ios /= 0) call gc_abort(me,NPROC,"ios error: forest.pcnt")
    endif

!hf u2    call stop_test(.true.,me,NPROC,ios,"ios error: forest.pcnt")

    if (me ==  0) then
!  
!      Read in 150km forest data (In percent of grid square),
!      normalise land-coverage so that sea areas are
!      corrected for. (The input file forest.pcnt only adds
!      up to one if there are no sea areas)
!
        do while (.true.)
          read(IO_FORES,*,err=1000,end=1000) itmp,jtmp,(tmp(i),i=1,6)
          itmp = 3*itmp + 34     ! Convert itmp from 150 to 50km grid
          jtmp = 3*jtmp + 10     ! Convert jtmp from 150 to 50km grid

          tmp(:) = 0.01 * tmp(:)
          sumland = sum( tmp )
          if( sumland >  1.0 ) tmp(:) = tmp(:)/sumland
          if ( DEBUG .and. any( tmp<0 )  ) then
             print *, "Biogenics: negative data at ", itmp, jtmp
	     inp_error = 999
             call gc_abort(me,NPROC,"ERROR:Forest_Init")
!hf u2	     goto 99    ! call stop_all("ERROR:Forest_Init")
          end if 

!          write(6,*) 'skogting',itmp,jtmp,sumland,(tmp(i),i=1,2)
!          Distribute percentage cover over 9 50km grids:
          do j = -1,1
              j50 = jtmp + j - JSMBEG + 1
              if(j50.ge.1.and.j50.le.GJMAX)then
           do i = -1,1
              i50 = itmp + i - ISMBEG + 1
              if(i50.ge.1.and.i50.le.GIMAX)then
!ds2.... add Grid area:   PLEASE CHECK MAP FACTOR USAGE !
!su            gridarea_m2 = 50000. * 50000. * xmd(i50,j50)
                    gforest(:,i50,j50) = tmp(:)
!                 do cat = 1, 6
!!                   gforest(cat,i50,j50) =  gridarea_m2 * tmp(cat)*invsland
!                    gforest(cat,i50,j50) = tmp(cat)
!                 enddo    !ncat
               endif        !i50
            enddo        !i
          endif        !j50
        enddo        !j
      enddo            !do while
1000  close(IO_FORES)
    endif  ! me = 0


!hf u299	call stop_test(.true.,me,NPROC,inp_error,"ERROR:Forest_Init")

!
!ds.......................end of new 1.................................
!
    call global2local(gforest,forest,MSG_READ1,6,GIMAX,GJMAX,1,1,1)

    forsum(:) = 0.0


    do j=lj0, lj1                !gv2 1,ljmax
       do i=li0, li1             !gv2 1,limax
            forsum(:)=forsum(:)+forest(:,i,j)
       end do
    end do
    call gc_rsum(6,NPROC,info,forsum)
    if (me == 0) write(unit=6,fmt="(a12,6f10.2)") "forest sum:", &
                                                  (forsum(cat),cat=1,6)



!ds  Calculate isoprene and terpene emitting biomass rates 
!     (at full sunlight, 30 deg.C)
!ds  e.g. oak   =  300.  *  41.2   *  forest(1,i,j)
!ds                g/m2  *  ug/g/h *    m2
!ds                =>  ug/m2/h

    isopsum = 0.0
    terpsum = 0.0

    do j = lj0, lj1     !gv2  1, ljmax   !gv  lj0,lj1
      do i = li0,li1    !gv2  1, limax   !gv  li0,li1

        !1) Isoprene
        oak   =  300.*41.2 *forest(1,i,j)
        decid =  300.* 1.2 *forest(2,i,j)
!ds         if(gb(i,j).lt.40.) then
!ds           decid = decid + 0.9*oak
!ds           oak = oak*0.1
!ds         end if
        conif = 1000.* 0.4 *forest(3,i,j)
        spruc = 1400.* 1.24*forest(4,i,j)
        crops =   40.* 0.3 *forest(5,i,j)
        sitka = 1400.* 6.24*forest(6,i,j)
!ds        Convert percent to fraction
!ds        xiso(i,j) = 0.01*(oak+decid+conif+spruc+crops+sitka)

       ! xiso -> emforest(BIO_ISOP
        emforest(BIO_ISOP,i,j) = (oak+decid+conif+spruc+crops+sitka)*xm2(i,j)

        !2) Terpene  (crude for Mexico job......)
        conif = 1000.* 3.0 *forest(3,i,j)
        spruc = 1400.* 3.0 *forest(4,i,j)
        sitka = 1400.* 3.0 *forest(6,i,j)

        emforest(BIO_TERP,i,j) = (conif+spruc+sitka)*xm2(i,j)

!       Exclude emissions from sea grids:  crude !!!
        if(iclass(i,j) == 0) emforest(:,i,j) = 0.0

        isopsum   = isopsum   + emforest(BIO_ISOP,i,j)
        terpsum   = terpsum   + emforest(BIO_TERP,i,j)

      end do
    end do

    call gc_rsum(1,NPROC,info,isopsum)
    call gc_rsum(1,NPROC,info,terpsum)
    if (me .eq. 0) then
      write(6,*) 'isopsum, terpsum',isopsum, terpsum
    end if

    !/-- Calculate canopy environmental correction factors
    !-----------------------------------------------------
    ! ag suffix  = Alex Guenther's parameter values

    agts = 303.
    agr = 8.314
    agtm = 314.
    agct1 = 95000.
    agct2 = 230000.
!
!ds     We have already isoprene emissions in ug/m2/h
!       6.023e23 - avogadros number (molecules/mol
!       1.e-6 --  micro-g ---> g 
!       1.e-4 --  m2 -> cm2
!       3600  -- hour-1 ----> s-1
!       68.  -- atw of isoprenes
!       ===>  ecf_isop no dimension * micro-g/m2/hour
!                                  to molecules/cm3/s
!
!    we multiply by chefac in setemis_ds
!
    fac = 1.e-9/(68.*3600.)
!    fac = 6.023e23*1.e-12/(68.*3600.)
!jej        fac = 6.023e23*1.e-10/(68.*3600.)       ! CHECK !!!!
!
!
!ds JOFFEN  I didn't understand the following line????
!ds     1.e-6 -- m2 ---> cm2 + m ---> cm (devided by box height later)
!ds     NB !!!!     product = 1
!
!ds     A simple assumption here that the light correction is approx. 1
!ds     (only applied in daylight I hope!!?)
!

    agcl = 1.0 ! this should be applied during daylight.

    if(me.eq.0) write(6,*) 'agcl',agcl

    do it = 1,40
      itk = it + 273.15
      agct = exp(agct1*(itk - agts)/(agr*agts*itk)) / &
                (1. + exp(agct2*(itk - agtm)/(agr*agts*itk)))

      canopy_ecf(BIO_ISOP,it) = agct*agcl*fac

      ! Terpenes
      agct = exp( 0.09*(itk-agts) )
      ! for terpene fac = 0.5*fac(iso): as mass terpene = 2xmass isoprene

      canopy_ecf(BIO_TERP,it) = agct * fac * 0.5   

      !if(me.eq.0) write(6,*) 'agcl',agcl,'agct',agct,ecf_isop(it)
      if(DEBUG .and. me.eq.0) &
             write(6,"(A12,i4,5g12.3)") 'Biogenic ecfs: ', &
                  it, agct, agcl, fac, &
                  canopy_ecf(BIO_ISOP,it), canopy_ecf(BIO_TERP,it)
    end do

   end subroutine Forests_Init
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module Biogenics_ml
