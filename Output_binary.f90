
module Output_binary_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     
!     Content:
!     
!     Output_binary - general output routine which sets the ident fields
!     for binary data, and calls the relevant 2D and 3D output subroutines:
!
!     Output2D  -  output of derived fields such as average air concentrations,
!     accumulated wet and dry deposition, etc.
!     
!     Output3D  -  output of derived fields in sigma levels
!
!    all the operations are now done in one call:
!        outarr_int2    cf. putflti2.F
!    here are only the calls and some presettings:
!        the specifications
!        the independent (characteristic for a given array)
!            scaling factor scale (to come i.e. to PPB)
!        total scale is then scal*10**(-iscal), iscal
!        is defined according to the maximum of the array
!        in outarr_int2 and written to ident(20), which is part
!        of the output fields
!        if the output array is some combination of other arrays,
!        then here this is assigned to rtmp    
!     
!----------------------------------------------------------------------------
!  Re-coded for unified model and Derived fields, Dave Simpson   16-17/10/01
!  Based upon MACHO version from jej.
!  Originally outchem, from:
!     Erik Berge, DNMI    Roar Skaalin, SINTEF Applied Mathematics
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Documentation f od DNMI/xfelt ident fields from jej, 5/01. 
! Code changed for position at the pole (xp, yp, read from infield).
! This is to make sure that the gridded data are located to the correct
! geographical location.
! ident(1) = producer  nr.   07 Washington, 88 Oslo, 98 ECMWF  ........
! ident(2) = grid area  (1841 for PARLAM PS)
! ident(3) = Data type
! ident(4) = time parameter (6 means a 6 hour forecast.... ) For most purposes
!   it does not matter what we use here. But for vertical slices it should be 
!   the same as the met. data (9 or 12). I dont know......
! ident(5) = vertical coordinat
!             = 1  P  (pressure coordinates)
!             = 2  sigma
!             = 3  h
!             = 4  theta
!             = 5 z (sea depth)
!             = 10  eta (hybrid)  used in standard HIRLAM
! ident(6) = parameter ???
! ident(7) = level 1 -- 20 or 1000 for surface
! ident(8) = level 2????    I have now set this to 0 (as it seemed to be in the
!    old days and in Krystoft made version. The value 980 seems to cause trouble,
!    and I dont understand where it comes from.
! ident(9) = grid type???
! ident(10) = # points in x
! ident(11) = # points in y
! ident(12) = year
! ident(13) = month + day    (ds- 100*nmonth+nday?)
! ident(14) = time*100+(min=0)
! ident(15) = x position of N pole
! ident(16) = y pole
! ident(17) = grid distance B*100
! ident(18) = l*100
! ident(19) + different things ?????
! ident(20) = scaling
!
 use Par_ml, only: GIMAX, GJMAX,       & ! x, y dimensions, (global domain)
                   me,NPROC            & ! processor number
                  ,MAXLIMAX, MAXLJMAX  & ! x,y dims on local domain
                  ,ISMBEG, JSMBEG        ! for restricted run-domain

 use ModelConstants_ml, only: identi, KMAX_MID

 use GridValues_ml, only: sigma_mid, xp, yp

 use Io_ml,           only : IO_OUT

 use My_Derived_ml, only : &
            NDDEP, NWDEP, NDERIV_2D, NDERIV_3D &
           ,Deriv                              & ! Derived data "type" 
           ,LENOUT2D, LENOUT3D                 & ! Dimension of averaging arrays
           ,IOU_INST, IOU_YEAR                 & ! Defines period of avg.
           ,IOU_MON, IOU_DAY                   & !  "   "
           ,f_wdep, f_ddep, f_2d, f_3d         & ! Definitions
           ,wdep, ddep, d_2d, d_3d               ! Data fields

 use Derived_ml, only : &
           nav_wdep, nav_ddep, nav_2d, nav_3d   ! No. terms for averging

 use NetCDF_ml,  only: Out_netCDF

 implicit none
 private


   public  :: Output_binary
   private :: Output2D
   private :: Output3D

  !/-- some variables 

   integer, private                ::  msgnr
   integer, dimension(20), private ::  ident

contains

  !<==========================================================================
   subroutine Output_binary(iotyp,outfilename)

     !   Output_binary - general output routine which sets the ident fields
     !   for binary data, and calls the relevant 2D and 3D output subroutines:
     !   
     !   the output parameters are as defined in My_Derived_ml,
     !   (old xddep,xwdep,xnav_m,xnshav_m,prodo3
     !      ,prdep,xnsu_m,aot40,aot60,accsu_m)
     !   
     !----------------------------------------------------

    integer ,      intent(in) :: iotyp
    character*30 , intent(in) :: outfilename

    integer i, k, n, ios

    ios = 0

    if(me == 0)then

      open(IO_OUT,file=outfilename,access='sequential',  &
              form='unformatted',status='unknown',iostat=ios)

      if(ios /= 0)then
        write(6,*) 'can not open outfile ',IO_OUT, 'filename : ',outfilename
        call gc_abort(me,NPROC,"error in outchem")
      endif

    endif

    ident(:) = identi(:)

    ident(4) = 6
    ident(5) = 2   ! ??? is this right ??
    ident(8) = 0   ! Set after Joffen's comments above
    ident(10) = GIMAX
    ident(11) = GJMAX

    if(identi(17) > 0)then !pw test for standard felt-format
       ident(15) = 100 * ( nint(xp)  - (ISMBEG - 1) )
       ident(16) = 100 * ( nint(yp)  - (JSMBEG - 1) )
    else
       ident(15) = ( nint(xp)  - (ISMBEG - 1) )  
       ident(16) = ( nint(yp)  - (JSMBEG - 1) )  
    endif


    if(me == 0) write(6,"(a12,2f6.1,/(10i6))") 'identi node 0', &
         xp, yp, (ident(k),k=1,20)



  ! Start output of data:
  !     
  !   put the 2-d derived fields, e.g. aot, accsu, So2, SO4, etc.
  !cccccccccccccccccccccccccc
     if ( NDERIV_2D > 0 ) &    ! u.1
      call Output2d(iotyp,NDERIV_2D, nav_2d, f_2d, d_2d)

  !  
  !     
  !***  Dry depositions
  !cccccccccccccccccccccccccc

     if ( NDDEP > 0 ) &    ! u.1
      call Output2d(iotyp,NDDEP, nav_ddep, f_ddep, ddep)


  !***  Wet deposition
  !cccccccccccccccccccccccccccc
  
     if ( NWDEP > 0 ) &    ! u.1
      call Output2d(iotyp,NWDEP, nav_wdep, f_wdep, wdep)


  !   put the 3-d derived fields, e.g. O3
  !cccccccccccccccccccccccccccc

     if ( NDERIV_3D > 0 ) &    ! u.1
      call Output3d(iotyp,NDERIV_3D, nav_3d, f_3d, d_3d)


      if (me == 0) close(IO_OUT)

   end subroutine Output_binary
!<==========================================================================

   subroutine  Output2d(iotyp, dim, nav, def, dat)

     integer,                         intent(in) :: iotyp
     integer,                         intent(in) :: dim ! No. fields
     integer, dimension(dim,LENOUT2D),intent(in) :: nav ! No. items averaged
     type(Deriv), dimension(dim),     intent(in) :: def ! definition of fields
     real, dimension(dim,MAXLIMAX,MAXLJMAX,LENOUT2D), &
                                      intent(in) :: dat ! Data arrays

     logical :: wanted     ! Set true for required year, month, day or inst.
     integer :: icmp       ! component index
     real    :: scale      ! Scaling factor
     character(len=size(f_2d%class)) :: typ  !  See defs of f_2d
     character*30 ::varname
     integer ::i,j,kmax ,ndim

       ident(7)  =  1000 ! Set artificial surface value for all 2-D fields


       do icmp = 1, dim
          typ = f_2d(icmp)%class
          wanted = .false.
          if( iotyp == IOU_YEAR) wanted = def(icmp)%year
          if( iotyp == IOU_MON ) wanted = def(icmp)%month
          if( iotyp == IOU_DAY ) wanted = def(icmp)%day
          if( iotyp == IOU_INST) wanted = def(icmp)%inst

          if ( wanted ) then 
              ident(6)  = def(icmp)%code      ! Identifier for DNMI/xfelt
              scale  = def(icmp)%scale
              if ( nav(icmp,iotyp) > 0  ) scale  = &
                                         def(icmp)%scale / nav(icmp,iotyp)
              msgnr = 5040 + icmp

              call outarr_int2(msgnr,IO_OUT,ident,dim,icmp, &
                                                 dat(1,1,1,iotyp),scale)
              if( iotyp == IOU_YEAR)then 
                 !writeout in netCDF output
                 write(varname,fmt='(''Out_2D'',i4.4)')def(icmp)%code
                 ndim=2                 
                 kmax=1
                 call Out_netCDF(ndim,varname,ident &
                   ,dim,kmax,icmp,dat(:,:,:,IOU_YEAR),scale)

              endif        

          endif        ! wanted
       enddo        !icmp

    end subroutine  Output2d
  !<==========================================================================

    subroutine  Output3d(iotyp, dim, nav, def, dat)

    ! Similar to 2-D data, but for now we allow an assumed-shape for the
    ! last dimension. This is because for 3-D ararys we currently only use 
    ! dimension IOU_INST:IOU_MON, instead of IOU_INST:IOU_DAY, to save
    ! memory space when using these big 3-D arrays.

     integer,                         intent(in) :: iotyp
     integer,                         intent(in) :: dim ! No. fields
     integer, dimension(dim,LENOUT3D), intent(in) :: nav ! No. items averaged
     type(Deriv), dimension(dim),     intent(in) :: def ! definition of fields
    !!   dimension(dim,MAXLIMAX,MAXLJMAX,KMAX_MID,LENOUT3D)
     real, dimension(:,:,:,:,:), &
                                      intent(in) :: dat ! Data arrays

     logical :: wanted     ! Set true for required year, month, day or inst.
     integer :: icmp, k    ! component index, vertical coord.
     real    :: scale      ! Scaling factor
     character*30 ::varname
     integer ::i,j,kmax,ndim 

       do icmp = 1, dim

          wanted = .false.
          if( iotyp == IOU_YEAR) wanted = def(icmp)%year
          if( iotyp == IOU_MON ) wanted = def(icmp)%month
          if( iotyp == IOU_DAY ) wanted = def(icmp)%day  ! Should be falseF
          if( iotyp == IOU_INST) wanted = def(icmp)%inst

          if ( wanted ) then 
              ident(6)  = def(icmp)%code      ! Identifier for DNMI/xfelt
              scale     = def(icmp)%scale


              if (iotyp /= IOU_INST  .and. & 
                   nav(icmp,iotyp) > 0  )  scale = def(icmp)%scale &
                                                 / nav(icmp,iotyp)


              do k=1,KMAX_MID
                 ident(7)  = k
                 ident(19) = sigma_mid(k)*10000.0 + 0.5

                 msgnr = 100*icmp + k

                 call outarr_int2(msgnr,IO_OUT,ident &
                   ,dim,icmp,dat(1,1,1,k,iotyp),scale)

              enddo    !k

              if( iotyp == IOU_YEAR)then 
                 !writeout in netCDF output
                 write(varname,fmt='(''Out_3D'',i4.4)')def(icmp)%code
                 ndim=3
                 kmax=KMAX_MID
                 call Out_netCDF(ndim,varname,ident &
                   ,dim,kmax,icmp,dat(:,:,:,:,IOU_YEAR),scale)

              endif        

            endif        ! wanted
          enddo        !icmp

    end subroutine  Output3d


 end module Output_binary_ml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
