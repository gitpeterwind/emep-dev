!TO DO:
!gc_send only before writing to disc
!interpolate bc (in time)
!more structured!
!could save only boundaries which really are outer boundaries
!write dates on disc?

   module Nest_ml
!
!This module performs the reading or writing of data for nested runs
!
  !
  !
  !1)Reads concentrations from disc. 
  !The concentrations stored on disc are defined in another 
  !coordinates system ("global" grid) and covers a larger area.
  !2)Convert coordinates
  !3)Interpolate the concentrations vertically and horizontally 
  !  to fit into the grid
  !
  !
  !
  !
  !
    use ModelConstants_ml,    only : KMAX_MID   ! vertical extent
    use GenSpec_adv_ml,  only: NSPEC_ADV         ! => No. species 
!    use GenSpec_shl_ml,  only: NSPEC_SHL         ! => No. species 
    use Par_ml   ,      only : MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,ISMBEG,JSMBEG &
                               , me, NPROC,li0,li1,lj0,lj1,limax,ljmax&
                               , tgi0, tgj0, tlimax, tljmax

  implicit none

  integer,parameter ::MODE=0   !0=donothing , 1=write , 2=read

  private

  !/-- subroutines

  public  :: readxn
  public  :: wrtxn

  private  :: readxnfromdisc
  private  :: updatebc
  private  :: GetGlobalData_Nest         ! Opens, reads bc_data, closes global data
  private  :: Getdataparameters
  private  :: InterpolationFactors_Nest  ! Gets factors for interpolation


  private :: vert_interpolation_Nest
  private :: lb2ij
  private :: ij2lb

  integer, parameter, private :: NGRIDSPECI=8
  real, private     ::gridspecifications(NGRIDSPECI)
  real, save, private::xpg,ypg,fig,ang,gridwidth_mg
  integer, save, private::kmaxg, gimaxg,gjmaxg

  integer, parameter :: NDATA=8     !how many data sets in each file
  integer, parameter :: NHOURSAVE=3 !time between two saves. should be a fraction of 24
  integer, parameter :: NHOURREAD=3 !time between two reads. should be a fraction of 24
                                    !and NHOURREAD >= NHOURSAVE

    real,save, dimension(NSPEC_ADV,MAXLIMAX,KMAX_MID,0:NDATA) :: xn_adv_bnds,xn_adv_bndn !north and south
    real,save, dimension(NSPEC_ADV,MAXLJMAX,KMAX_MID,0:NDATA) :: xn_adv_bndw,xn_adv_bnde !west and east

contains

  subroutine readxn(indate)
    use Dates_ml,       only : date     ! No. days per year, date-type 
    type(date), intent(in) :: indate           ! Gives year..seconds
    integer,save  :: NDATAcount=0,first_data=0
    character*30  :: filename


    integer :: errcode,io_num,datafound

    if(MODE /= 2)return
     if(me==0)   print *,'call to READXN',me
    if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
    if(first_data==-1)then
       first_data=0
       return
    endif
    if(mod(indate%hour,NHOURSAVE)/=0)return
    NDATAcount=NDATAcount+1
    if(me==0)print *,'READXN:  NDATAcount=',NDATAcount
    if(NDATAcount==1)then 
       !fetch a new dataset from disc
       write(filename,fmt='(''NestDATA_50km'',3i2.2)')indate%month,indate%day,indate%hour
       if(me==0)print *,'READXN: will try to find data in ',filename
       call readxnfromdisc(filename,datafound)
       if(datafound==0)then
          !no data was found. continue program without updating BC
          if(me==0)then
          print *,'WARNING: READXN: no datafile was found ',filename
          endif
          NDATAcount=0
          return
       endif
    endif
    
    if(mod(indate%hour,NHOURREAD)==0)then
    if(me==0)print *,'READXN: update boundaries ',NDATAcount,indate%hour
    call updatebc(NDATAcount) !update boundary concentrations
    endif

    if(NDATAcount==NDATA)NDATAcount=0
    
    return
  end subroutine readxn

  subroutine wrtxn(indate)

    !to be  "serialized"

    use GenSpec_adv_ml, only : NSPEC_ADV         ! => No. species 
    use Chemfields_ml,  only : xn_adv, xn_shl    ! emep model concs.
    use Io_ml   ,       only : IO_NEST
    use GridValues_ml,  only : xp,yp,fi,an,GRIDWIDTH_M
    use Dates_ml,       only : date     ! No. days per year, date-type 

    implicit none

    type(date), intent(in) :: indate           ! Gives year..seconds


    integer,save  :: NDATAcount=0
    integer  :: n,i,j,k,i0,j0,i1,j1,i2,j2,d,imaxs,jmaxs,alloc_err
    integer  :: iminsgrid,imaxsgrid,jminsgrid,jmaxsgrid
    integer  :: info,msnr1,msnr2
    real, allocatable,dimension(:,:)  :: sl,sb,ir,jr
    real, allocatable,dimension(:,:,:,:)  :: xn_adv_temp
    real, save, allocatable,dimension(:,:,:,:,:)  :: xn_adv_store
    real  :: fis,ans,xps,yps
    character*30, save  :: filename
    integer, save  :: my_first_call=0

    if(MODE /= 1)return

!    print *,'Nest_ml',me,indate%day,indate%hour,indate%seconds
    if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
!    if(mod(indate%hour,NHOURSAVE)/=0)return

!    print *,'save data',me,indate%day,indate%hour


    ! small grid parameters:

    imaxs = 75
    jmaxs = 60
    fis   = 10.50000
    ans   = 11888.44824218750
    xps   = 41.006530761718750
    yps   = 3234.5815429687500

    imaxs = 7
    jmaxs = 7
    fis   = -32
    ans   = 237.73164421375
    xps   = 43.-82.5
    yps   = 121.-70.5

    imaxs = 110
    jmaxs = 110
    fis   = -32
    ans   = 5*237.73164421375
    xps   = -165.+2.
    yps   = 285.+2.

!    imaxs = 115
!    jmaxs = 115
!    fis   = -32
!    ans   = 10*237.73164421375
!    xps   = -519.5+1.
!    yps   = 620.5+1.

    allocate(sl(imaxs,jmaxs), stat=alloc_err)
    allocate(sb(imaxs,jmaxs), stat=alloc_err)
    allocate(ir(imaxs,jmaxs), stat=alloc_err)
    allocate(jr(imaxs,jmaxs), stat=alloc_err)
    if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "alloc failed in writeconcentrations")


    ! find the area that needs to be stored
    ! 1) find longitudes and latitudes of small s (local) grid

    call ij2lb(imaxs,jmaxs,sl,sb,fis,ans,xps,yps)

    ! 2) find the corresponding i,j in large ("global g") map 
    !    coordinates will be real numbers ir,jr

    call lb2ij(imaxs,jmaxs,sl,sb,ir,jr,fi,an,xp,yp)


    ! 3) find a rectangle in the large grid that contains all the coordinates of the small grid

    iminsgrid = int(minval(ir(1:imaxs,1:jmaxs)))
    imaxsgrid = int(maxval(ir(1:imaxs,1:jmaxs)))+1
    jminsgrid = int(minval(jr(1:imaxs,1:jmaxs)))
    jmaxsgrid = int(maxval(jr(1:imaxs,1:jmaxs)))+1

    if(my_first_call==0.and.me==0)then
    write(*,*)minval(ir(1:imaxs,1:jmaxs))
    write(*,*)minval(ir(1:imaxs,1:jmaxs))
    write(*,*)int(minval(ir(1:imaxs,1:jmaxs))-0.00001)
    write(*,*)'small grid ',minval(ir(1:imaxs,1:jmaxs)),minval(jr(1:imaxs,1:jmaxs))
    write(*,*)'small grid ',iminsgrid,imaxsgrid,jminsgrid,jmaxsgrid
    write(*,*)'large grid ',ISMBEG,ISMBEG+GIMAX-1,JSMBEG,JSMBEG+GJMAX-1
    endif

    deallocate(sl, stat=alloc_err)
    deallocate(sb, stat=alloc_err)
    deallocate(ir, stat=alloc_err)
    deallocate(jr, stat=alloc_err)
    if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in writeconcentrations")

    ! coordinates in the domain
    i1 = iminsgrid-ISMBEG+1
    !     i2 = imaxsgrid-ISMBEG+1
    i2 = imaxsgrid-ISMBEG+1-i1+1 !number of cells in i direction
    j1 = jminsgrid-JSMBEG+1
    !     j2 = jmaxsgrid-JSMBEG+1
    j2 = jmaxsgrid-JSMBEG+1-j1+1 !number of cells in j direction

    if(my_first_call==0.and.me==0)then
       write(*,*)'i,j of small grid ',i1,i2,j1,j2
       my_first_call=1
    endif

    ! check if the small grid is inside the large grid
    if(i1<1 .or. i1+i2-1 > GIMAX .or. j1<1 .or. j1+j2-1 > GJMAX )then
       write(*,*)'The small grid is not inside the large grid!'
       write(*,*)'i,j of small grid ',i1,i2,j1,j2
       write(*,*)iminsgrid,imaxsgrid,jminsgrid,jmaxsgrid
       write(*,*)ISMBEG,ISMBEG+GIMAX-1,JSMBEG,JSMBEG+GJMAX-1
       call gc_abort(me,NPROC, "cannot write down required concentrations")
    endif

    ! 4) gather the data from other processors
    if(me.eq.0)then
       allocate(xn_adv_temp(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID), stat=alloc_err)
       if(NDATAcount==0)then
          allocate(xn_adv_store(NSPEC_ADV,i2,j2,KMAX_MID,NDATA), stat=alloc_err)
       endif
    if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "alloc failed in writeconcentrations 4")

    endif

    if(me==0.and.NDATAcount==0)then 
       write(filename,fmt='(''NestDATA_50km'',3i2.2)')indate%month,indate%day,indate%hour
       print *,'DEFINING FILE ',filename
    endif
    if(me.eq.0)NDATAcount=NDATAcount+1
       
    do d = 1,NPROC-1
!       if(me.eq.0)             write(*,*)d,tgi0(d),tlimax(d),tgj0(d),tljmax(d)

       if(tgi0(d) <= i1+i2-1 .and. tgi0(d) + tlimax(d)-1 >= i1 .and.&
            tgj0(d) <= j1+j2-1 .and. tgj0(d) + tljmax(d)-1 >= j1) then

          ! subdomain contains a part of the small grid 
    msnr1 = 88+d
    msnr2 = 88+NPROC+d

          if(me == d )then
             call gc_rsend(msnr1, NSPEC_ADV*MAXLIMAX*MAXLJMAX*KMAX_MID &
                           , 0, info, xn_adv_temp, xn_adv)	  
          endif
          if(me == 0 )then
             call gc_rrecv(msnr1,NSPEC_ADV*MAXLIMAX*MAXLJMAX*KMAX_MID &
                              , d, info, xn_adv_temp, xn_adv)
             do j = 1, tljmax(d)
                j0 = tgj0(d)+j-j1
                if( j0 >= 1 .and. j0 <= j2 ) then
                   do i = 1,tlimax(d)
                      i0 = tgi0(d)+i-i1
                      if( i0 >= 1 .and. i0 <= i2 ) then
!                         write(*,*)d,i,j,i0,j0,i+tgi0(d)-1,j+tgj0(d)-1
                         do k=1,KMAX_MID
                            do n=1,NSPEC_ADV
                               xn_adv_store(n,i0,j0,k,NDATAcount)=xn_adv_temp(n,i,j,k)
                            enddo
                         enddo
                      endif
                   enddo
                endif
             enddo
          endif
       endif
    enddo

    if(me == 0 )then
    deallocate(xn_adv_temp, stat=alloc_err)

    if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in writeconcentrations")
    endif

! also include xn from me=0 

    if(me == 0 )then
    do d = 0,0

       if(tgi0(d) <= i1+i2-1 .and. tgi0(d) + tlimax(d)-1 >= i1 .and.&
            tgj0(d) <= j1+j2-1 .and. tgj0(d) + tljmax(d)-1 >= j1) then
          
          ! subdomain contains a part of the small grid 
          
          do j = 1, tljmax(d)
             j0 = tgj0(d)+j-j1
             if( j0 >= 1 .and. j0 <= j2 ) then
                do i = 1,tlimax(d)
                   i0 = tgi0(d)+i-i1
                   if( i0 >= 1 .and. i0 <= i2 ) then
                      write(*,*)d,i,j,i0,j0
                      do k=1,KMAX_MID
                         do n=1,NSPEC_ADV
                            xn_adv_store(n,i0,j0,k,NDATAcount)=xn_adv(n,i,j,k)
                         enddo
                      enddo
                   endif
                enddo
             endif
          enddo
       endif
    enddo
    endif

    ! 5) write the relevant part of the concentrations onto disc

    if(me.eq.0)then
!       print *,xn_adv_store(2,i0,j0,17,NDATAcount)
       if(NDATAcount==NDATA)then

       write(*,*)'write new NEST file: ',filename
       write(*,*)'dimensions: ',NSPEC_ADV, i2,j2,KMAX_MID
       write(*,*)'coordinates: ',xp-ISMBEG-i1+2,yp-JSMBEG-j1+2,fi,an,GRIDWIDTH_M
       open(unit=IO_NEST,file=filename,form='unformatted',action='write')

       write(IO_NEST)NSPEC_ADV, i2,j2,KMAX_MID
       write(IO_NEST)xp-ISMBEG-i1+2,yp-JSMBEG-j1+2,fi,an,GRIDWIDTH_M

       write(IO_NEST)xn_adv_store

       close(IO_NEST)

       if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in writeconcentrations")

       NDATAcount=0
       deallocate(xn_adv_store, stat=alloc_err)
       endif
    endif
    return
  end subroutine wrtxn


  subroutine readxnfromdisc(filename,datafound)

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !** DESCRIPTION
    !   read in monthly-average global mixing ratios, and if found, collect the 
    !   data in  bc_adv for later interpolations
    !   (NB!!  if mixing ratio by mass the scale by molcular weight)
    !   ds- comment - so far no scaling is done, but this could be done
    !   in Set_bcmap with atomic weights.... for the future..
    !
    ! On the first call, we also run the setup-subroutines
    !____________________________________________________________________________
    use Chemfields_ml,         only: xn_adv  ! emep model concs.
    use GridValues_ml,         only: gl, gb    ! lat, long

    implicit none

    integer, intent(inout) :: datafound
    character*30, intent(in)  :: filename

    integer :: k,i,j     ! loop variables
    integer :: info              !  used in gc_rsend
    integer :: io_num            !  i/o number used for reading global data

    integer alloc_err1,alloc_err2,alloc_err3,alloc_err,NDATAcount

    real, allocatable,dimension(:,:,:,:,:) :: bc_adv

    real,save, dimension(MAXLIMAX,MAXLJMAX) :: wt_00, wt_01, wt_10, wt_11
    integer,save, dimension(MAXLIMAX,MAXLJMAX) :: &
         ixp, iyp !  global model coordinates of point (i,j)


    integer  ::  errcode,n
    integer,save  ::reset3d=0
    


    if(me == 0)then
!       print *, "READXN: find dataparameters"
       call Getdataparameters(filename,errcode)

       gridspecifications = &
            (/ xpg,ypg,fig,ang,gridwidth_mg,real(kmaxg),real(gimaxg),real(gjmaxg) /)

!       write(*,*)'gridspecifications ',me,xpg,ypg,fig,ang,gridwidth_mg,kmaxg,gimaxg,gjmaxg
!       write(*,*)'dimensions ',NSPEC_ADV ,gimaxg,gjmaxg,KMAX_MID

    endif

    call gc_rbcast(416, NGRIDSPECI,0,NPROC,info,gridspecifications)

    xpg = gridspecifications(1)
    ypg = gridspecifications(2)
    fig = gridspecifications(3)
    ang = gridspecifications(4)
    gridwidth_mg = gridspecifications(5)
    kmaxg = nint(gridspecifications(6))
    gimaxg = nint(gridspecifications(7))
    gjmaxg = nint(gridspecifications(8))

    if(gridwidth_mg==0)then
       !the file was not found
       datafound=0
       return
    endif
    datafound=1

    allocate(bc_adv(NSPEC_ADV ,gimaxg,gjmaxg,KMAX_MID,NDATA), stat=alloc_err)
    if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "bc_adv alloc failed")


    if(me == 0)then

       call GetGlobalData_Nest(bc_adv,errcode)
       if(errcode /= 0 ) call gc_abort(me,NPROC,"ERROR Nest: GetGlobalData")

    endif


    call gc_rbcast(417, NSPEC_ADV*gimaxg*gjmaxg*KMAX_MID*NDATA,0,  &
         NPROC,info,bc_adv)


    call InterpolationFactors_Nest(gb,gl,ixp,iyp,wt_00,wt_01,wt_10,wt_11)

! check if the corners of the domain are inside the area covered by the large grid:
! (In principle we should test for all i,j , but test the corners should be good 
!  enough in practice) 

     if(int(ixp(1,1)) < 1 .or. int(ixp(1,1))+1 > gimaxg .or. &
        int(ixp(limax,1)) < 1 .or. int(ixp(limax,1))+1 > gimaxg .or. &
        int(ixp(1,ljmax)) < 1 .or. int(ixp(1,ljmax))+1 > gimaxg .or. &
        int(ixp(limax,ljmax)) < 1 .or. int(ixp(limax,ljmax))+1 > gimaxg .or. &
        int(iyp(1,1)) < 1 .or. int(iyp(1,1))+1 > gjmaxg .or. &
        int(iyp(limax,1)) < 1 .or. int(iyp(limax,1))+1 > gjmaxg .or. &
        int(iyp(1,ljmax)) < 1 .or. int(iyp(1,ljmax))+1 > gjmaxg .or. &
        int(iyp(limax,ljmax)) < 1 .or. int(iyp(limax,ljmax))+1 > gjmaxg ) then
       write(*,*)'Did not find all the necessary concentrations in file'
       write(*,*)'values needed: '
       write(*,*)ixp(1,1),iyp(1,1),wt_00(1,1),wt_01(1,1),wt_10(1,1),wt_11(1,1)
       write(*,*)ixp(limax,1),iyp(limax,1)
       write(*,*)ixp(1,ljmax),iyp(1,ljmax)
       write(*,*)ixp(limax,ljmax),iyp(limax,ljmax)
       write(*,*)'max values found: ',gimaxg ,gjmaxg
       call gc_abort(me,NPROC, "Nest3d: area too small")
    endif

    !===================================
    NDATAcount=1
    if(reset3d==0)then
       !reset concentrations on the whole domain
       print *, "NEST: RESET 3D ", me
    do k = 1, KMAX_MID
       do j = 1, ljmax
          do i = 1,limax
             do n=1,NSPEC_ADV
                xn_adv(n,i,j,k) =  &
                     wt_00(i,j) * bc_adv(n,ixp(i,j), iyp(i,j), k,NDATAcount) +  & 
                     wt_01(i,j) * bc_adv(n,ixp(i,j), iyp(i,j)+1,k,NDATAcount) +  & 
                     wt_10(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j), k,NDATAcount) +  & 
                     wt_11(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j)+1,k,NDATAcount)
             enddo
          enddo
       enddo
    enddo
    reset3d=1
    endif

 !save concentrations at boundaries only
    !save last "old" concentrations
    xn_adv_bnds(:,:,:,0)=xn_adv_bnds(:,:,:,NDATA)
    xn_adv_bndn(:,:,:,0)=xn_adv_bndn(:,:,:,NDATA)
    xn_adv_bndw(:,:,:,0)=xn_adv_bndw(:,:,:,NDATA)
    xn_adv_bnde(:,:,:,0)=xn_adv_bnde(:,:,:,NDATA)
    !make and save new boundaries concentrations
    do k = 1, KMAX_MID
       do j = 1, lj0
          do i = 1,limax
             do n=1,NSPEC_ADV
                xn_adv_bnds(n,i,k,1:NDATA) =  &
                     wt_00(i,j) * bc_adv(n,ixp(i,j), iyp(i,j), k,:) +  & 
                     wt_01(i,j) * bc_adv(n,ixp(i,j), iyp(i,j)+1,k,:) +  & 
                     wt_10(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j), k,:) +  & 
                     wt_11(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j)+1,k,:)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = lj1, ljmax
          do i = 1,limax
             do n=1,NSPEC_ADV
                xn_adv_bndn(n,i,k,1:NDATA) =  &
                     wt_00(i,j) * bc_adv(n,ixp(i,j), iyp(i,j), k,:) +  & 
                     wt_01(i,j) * bc_adv(n,ixp(i,j), iyp(i,j)+1,k,:) +  & 
                     wt_10(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j), k,:) +  & 
                     wt_11(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j)+1,k,:)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = 1, ljmax
          do i = 1,li0
             do n=1,NSPEC_ADV
                xn_adv_bndw(n,j,k,1:NDATA) =  &
                     wt_00(i,j) * bc_adv(n,ixp(i,j), iyp(i,j), k,:) +  & 
                     wt_01(i,j) * bc_adv(n,ixp(i,j), iyp(i,j)+1,k,:) +  & 
                     wt_10(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j), k,:) +  & 
                     wt_11(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j)+1,k,:)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = 1, ljmax
          do i = li1,limax
             do n=1,NSPEC_ADV
                xn_adv_bnde(n,j,k,1:NDATA) =  &
                     wt_00(i,j) * bc_adv(n,ixp(i,j), iyp(i,j), k,:) +  & 
                     wt_01(i,j) * bc_adv(n,ixp(i,j), iyp(i,j)+1,k,:) +  & 
                     wt_10(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j), k,:) +  & 
                     wt_11(i,j) * bc_adv(n,ixp(i,j)+1,iyp(i,j)+1,k,:)
             enddo
          enddo
       enddo
    enddo


    !===================================


    deallocate(bc_adv,stat=alloc_err2)
    if ( alloc_err2 /= 0 ) call gc_abort(me,NPROC,"de-alloc_err2")



  end subroutine readxnfromdisc


  subroutine GetGlobalData_Nest(bc_adv,errcode)

    use Io_ml,           only: IO_NEST, ios, open_file

	implicit none


    real, dimension(NSPEC_ADV,gimaxg,gjmaxg,KMAX_MID,NDATA), &
         intent(out) :: bc_adv   
    integer,            intent(inout) :: errcode    !  i/o number

    real, allocatable,dimension(:,:,:,:,:) :: bc_rawdata_adv   

    character(len=30) :: fname    ! input filename
    integer, save     :: oldmonth = -1  
    logical, save                    :: my_first_call = .true.
    integer :: i,j,alloc_err=0
    integer :: nspec_advg


    errcode = 0

!       write(unit=fname,fmt="(a6,i2.2,a4)") "gl_ass",month,".dat"
!       fname='concfile'
!    open(unit=IO_NEST,file=fname,form='unformatted',action='read')

      allocate(bc_rawdata_adv(NSPEC_ADV ,gimaxg,gjmaxg,kmaxg,NDATA), &
         stat=alloc_err)
      if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "bc_raw_adv alloc failed")

    read(IO_NEST) bc_rawdata_adv
    call vert_interpolation_Nest(NSPEC_ADV,bc_rawdata_adv,bc_adv)

    close(IO_NEST)
    
    deallocate(bc_rawdata_adv,stat=alloc_err)
    if ( alloc_err /= 0 ) call gc_abort(me,NPROC,"adv de-alloc_err")


  end subroutine GetGlobalData_Nest


  subroutine Getdataparameters(fname,errcode)

    use Io_ml,           only: IO_NEST, ios, open_file

	implicit none

    integer :: io_num    !  i/o number
    integer,            intent(inout) :: errcode    !  i/o number

    character(len=30) :: fname    ! input filename

    logical, save                    :: my_first_call = .true.
    integer :: i,j,alloc_err=0
    integer :: nspec_advg

    io_num = IO_NEST              ! for closure in BoundCOnditions_ml

    errcode = 0

!       write(unit=fname,fmt="(a6,i2.2,a4)") "gl_ass",month,".dat"
!       fname='concfile'
       write(*,*)'opening ', fname
       open(unit=IO_NEST,file=fname,form='unformatted',action='read',status='old',err=90)
       write(*,*)'opening file succesful: ', fname

       read(IO_NEST,err=90)nspec_advg, gimaxg,gjmaxg,kmaxg

       read(IO_NEST,err=90)xpg,ypg,fig,ang,gridwidth_mg

       if(nspec_advg.ne.NSPEC_ADV)then
          write(*,*)'wrong number of species! ',nspec_advg
          write(*,*)'should be ',NSPEC_ADV
          call gc_abort(me,NPROC, "wrong nspec")
       endif
       if(kmaxg.ne.20)then
          write(*,*)'20 levels hardcoded ',kmaxg
          call gc_abort(me,NPROC, "wrong kmaxg")
       endif
       return
90 continue
       errcode = 1

       xpg=0
       ypg=0
       fig=0
       ang=0
       gridwidth_mg=0
      write(*,*)'file does not exist ',fname
  end subroutine Getdataparameters


  subroutine vert_interpolation_Nest(nspec,bc_rawdata,bc_data)

    use ModelConstants_ml, only: PT                         ! Top of EMEP model (=1.0e4 Pa)
    use GridValues_ml, only: sigma_mid    ! EMEP sigma values

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !+
    ! Calculates interpolation from global-model sigma to EMEP model sigma. If
    ! sigma in emep model is greater than sigma in lowest layer in global model,
    ! concentrations frm the lowest global layer are used. 
    !    This routine  is run every time new bc data are read
    ! in. We could save the interpolation factors, but the calculations are
    ! trivial so we don't bother.
    !
    !NB: We assume that PT is the same in the two definitions of sigma levels
    !
    !--------------------------------------------------------------------------
	implicit none

    integer , intent(in)  :: nspec
    real,  dimension(nspec,gimaxg,gjmaxg,kmaxg,NDATA), intent(in)  :: bc_rawdata !  data 
    real,  dimension(nspec,gimaxg,gjmaxg,KMAX_MID,NDATA), intent(out) :: bc_data    !  data 

    logical, save                    :: my_first_call = .true.
    integer, dimension(KMAX_MID), save  :: kglob1, kglob2    !  k-levels from global model
    real,    dimension(KMAX_MID), save  :: c   ! Interpolation factor


    real ::  midsiga     &  ! global sigma value centered imediately below..(sigma larger)
         , midsigb        ! global sigma value centered imediately above..(sigma smaller)

    real    ::  a, b       ! help variables
    integer ::  k, kg, ish ! help variables
    integer i,j,i1,j1,iglbegw,iglendw,jglbeg,jglend

    real, dimension(20) ::  mygsig,sigmamid20
    real :: sigmamid17(17)


    if ( my_first_call ) then
       !-consistency check
       if ( size(bc_rawdata,4) /= 20        ) then  
          call gc_abort(me,NPROC, "wrong kmaxg")
       end if

       sigmamid17 = (/  &
            0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.76, 0.8475, 0.8975,&
            0.932, 0.9545,  0.9725,  0.985,  0.993,  0.997,  0.999  /)
       sigmamid20 = (/  &
            0.02, 0.06, 0.1, 0.1425, 0.195, 0.2635, 0.347, 0.4365, 0.5215, 0.599,&
            0.6695, 0.733,  0.7895,  0.839,  0.8815,  0.917,  0.9455, 0.967, 0.982, 0.994  /)

!       kmaxg = 20

       mygsig(1:kmaxg) = sigmamid20(1:kmaxg)


       ! - set up vertical interpolation factors, c, plus coordinate indices 
       !      kglob,ish
       !pw    midsiga(kg) > sigma_mid(k) > midsigb(kg) (except possibly for the lowest levels)

       kg = kmaxg-1
       midsigb = mygsig(kg)
       midsiga = mygsig(kg+1)

       do k = KMAX_MID,1,-1    ! size = KMAX_MID

          do while (sigma_mid(k) < midsigb .and. kg>1 )
             kg = kg - 1
             midsiga = midsigb
             midsigb = mygsig(kg)
          enddo

          if(sigma_mid(k) > midsiga)then
             c(k) = 1.  
          else
             c(k) = (midsigb-sigma_mid(k))/(midsigb - midsiga)
          endif

          kglob1(k) = kg       ! Sets global k-coordinate corresponding to k

       end do ! k
    end if  ! my_first_call
    !=======================
    my_first_call = .false.
    !=======================

    iglbegw=1
    iglendw=gimaxg
    jglbeg=1
    jglend=gjmaxg

    do k = 1,KMAX_MID
       j1 = 1
       kg = kglob1(k)
!       write(*,*)k,kg,c(k),sigmamid17(k),sigmamid20(kg)
!       write(*,*)sigmamid17(k),c(k) * sigmamid20(kg+1)+(1.0-c(k))* sigmamid20(kg)
       do j = jglbeg,jglend
          i1 = 1
          do i = iglbegw,iglendw

             bc_data(:,i1,j1,k,:) =  c(k) * bc_rawdata(:,i,j,kg+1,:) +  & 
                  (1.0-c(k))* bc_rawdata(:,i,j,kg,:)
             i1 = i1+1
          enddo
          j1 = j1+1
       enddo
    end do ! k


  end subroutine vert_interpolation_Nest
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine InterpolationFactors_Nest(gb,gl,ixp,iyp,wt_00, wt_01, wt_10, wt_11)

    use GridValues_ml , only : GRIDWIDTH_M
    use Functions_ml, only: bilin_interpolate  ! bilinear interpolation

	implicit none

	real, dimension(MAXLIMAX, MAXLJMAX), intent(in)::gb,gl
	integer, dimension(MAXLIMAX, MAXLJMAX), intent(inout)::ixp,iyp
	real, dimension(MAXLIMAX, MAXLJMAX), intent(out)::wt_00, wt_01, wt_10, wt_11
        real, dimension(MAXLIMAX, MAXLJMAX) ::x,y
        real ::ismbegg,jsmbegg



    !x, y : coordinates in the large ("global") map
    !               fig = 10.50000
    !               ang = 11888.44824218750
    !               xpg = 41.006530761718750
    !               ypg = 3234.5815429687500
!    fig = -32.0                     ! read from file
!    ismbegg=80.
!    jsmbegg=70.
!        ang = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.   
!    ang = an*GRIDWIDTH_M /GRIDWIDTH_Mg   
!    xpg = 43. -ismbegg+1.                
!    ypg = 121. -jsmbegg+1.               

    call lb2ij(MAXLIMAX, MAXLJMAX, gl,gb,x,y,fig,ang,xpg,ypg)


    !============================================================
    call bilin_interpolate(x,y,ixp,iyp,wt_00,wt_01,wt_10,wt_11)
    !============================================================

  end subroutine InterpolationFactors_Nest




  subroutine lb2ij(imax,jmax,gl,gb,ir2,jr2,fi2,an2,xp2,yp2)
    !-------------------------------------------------------------------! 
    !      calculates coordinates ir2, jr2 (real values) from gl(lat),gb(long) 
    !
    !      input:  xp2,yp2:   coord. of the polar point in grid2
    !              an2:   number of grid-distances from pole to equator in grid2.
    !              fi2:      rotational angle for the grid2 (at i2=0).
    !              i1max,j1max: number of points (grid1) in  x- og y- direction
    !
    !
    !      output: i2(i1,j1): i coordinates in grid2 
    !              j2(i1,j1): j coordinates in grid2 
    !-------------------------------------------------------------------! 


    implicit none


    integer :: imax,jmax,i1, j1
    real    :: fi2,an2,xp2,yp2
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: ir2(imax,jmax),jr2(imax,jmax)

    real, parameter :: PI=3.14159265358979323
    real    :: PId4,dr,dr2


    PId4    = PI/4.      
    dr2    = PI/180.0/2.      ! degrees to radians /2
    dr    = PI/180.0      ! degrees to radians 

    do j1 = 1, jmax
       do i1 = 1, imax

          ir2(i1,j1)=xp2+an2*tan(PId4-gb(i1,j1)*dr2)*sin(dr*(gl(i1,j1)-fi2))
          jr2(i1,j1)=yp2-an2*tan(PId4-gb(i1,j1)*dr2)*cos(dr*(gl(i1,j1)-fi2))
!          write(*,*)i1,j1,ir2(i1,j1),jr2(i1,j1),gl(i1,j1),gb(i1,j1)

       end do ! i
    end do ! j

    return
  end subroutine lb2ij

  subroutine ij2lb(imax,jmax,gl,gb,fi,an,xp,yp)
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              an:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid (at i=0).
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): longitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): latitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 



    implicit none


    integer :: i, j, imax, jmax
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: fi, an, xp, yp 
    real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
    real, parameter :: PI=3.14159265358979323


!    fi = -32.0
    glmin = -180.0

    glmax = glmin + 360.0
    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, jmax          
       dy  = yp - j            
       dy2 = dy*dy
       do i = 1, imax       

         dx = i - xp    ! ds - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl(i,j)=rl                     !     longitude
         gb(i,j)=rb                     !     latitude
!         write(*,*)i,j,gl(i,j),gb(i,j)
       end do ! i
    end do ! j

   return
  end subroutine ij2lb

  subroutine updatebc(NDATAcount) !update boundary concentrations

    use Chemfields_ml,  only : xn_adv    ! emep model concs.

    implicit none

    integer :: i,j,k,n,NDATAcount

    do k = 1, KMAX_MID
       do j = 1, lj0
          do i = 1,limax
             do n=1,NSPEC_ADV
                xn_adv(n,i,j,k) =  xn_adv_bnds(n,i,k,NDATAcount)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = lj1, ljmax
          do i = 1,limax
             do n=1,NSPEC_ADV
                xn_adv(n,i,j,k) =  xn_adv_bndn(n,i,k,NDATAcount)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = 1, ljmax
          do i = 1,li0
             do n=1,NSPEC_ADV
                xn_adv(n,i,j,k) =  xn_adv_bndw(n,j,k,NDATAcount)
             enddo
          enddo
       enddo
    enddo
    do k = 1, KMAX_MID
       do j = 1, ljmax
          do i = li1,limax
             do n=1,NSPEC_ADV
                xn_adv(n,i,j,k) =  xn_adv_bnde(n,j,k,NDATAcount)
             enddo
          enddo
       enddo
    enddo
!    print *,'BC updated ',me
  end subroutine updatebc

end module Nest_ml

