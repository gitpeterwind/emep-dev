!
!{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
!
   Module ReadField_ml
!
!{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
!
! -----------------------------------------------------------
! Reads iascii  real and integer fields, usually for the whole 
! model area, and calls global2local to distribute these to
! the calling processor. Fields initialised
!
! Note that only one value per grid square is allowed for so
! far.
!
! Written October 2001, HF
! BIG introduced, ds, Nov. 2001
!------------------------------------------------------------
  use Par_ml,        only : IILARDOM,JJLARDOM      &
                      ,ISMBEG,JSMBEG               &
                      ,MAXLIMAX,MAXLJMAX           &
                      ,MSG_READ7,MSG_READ5         &
                      ,me,NPROC                    &
                      ,GIMAX,GJMAX
  use Io_ml,         only : ios, open_file
  use GridValues_ml, only : fi, an, xp, yp, ij2ij !pw emep1.2beta

  integer,private  :: i, j, n         ! Local variables
  integer, private, parameter :: BIG = IILARDOM*JJLARDOM


  interface ReadField
     module procedure ReadField_r, ReadField_i
  end interface

contains

!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
  subroutine ReadField_r(IO_INFILE,fname,local_field)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character*20, intent(in) :: fname         ! File name
  real, intent(out) :: local_field(MAXLIMAX,MAXLJMAX)! Local field

  !/-- We need an array for the whole model domain

  real :: in_field(IILARDOM,JJLARDOM)! Field to be read
  real :: out_field(GIMAX ,GJMAX )!pw

  integer :: imaxin,jmaxin,imaxout,jmaxout,i0,j0 !pw
  real :: fi_EMEP,an_EMEP,xp_EMEP,yp_EMEP,fiout,anout,xpout,ypout !pw

    if (me==0)then
       call open_file(IO_INFILE,"r",fname,needed=.true.)
          if ( ios /= 0 )then
          write(6,*) 'error in opening IO_INFILE', fname
          call gc_abort(me,NPROC,"newmonth: error opening in_field")
          endif
    endif !me==0


    if (me == 0) then

          READFIELD : do n = 1, BIG
             read(IO_INFILE,*,iostat=ios) i,j,in_field(i,j)
             if (  i < 1 .or. i > IILARDOM  .or. &
                   j < 1 .or. j > JJLARDOM  ) then  
                  write(6,*)'error in i,j index in IO_INFILE=',fname, i,j
                  ios = 88  !! Set an error code
                  exit READFIELD
             endif
            if ( ios /= 0 ) exit READFIELD
          enddo READFIELD

       close(IO_INFILE)
       if ( ios /= 0 )then             
           write(6,*) 'error in reading IO_INFILE', fname
           call gc_abort(me,NPROC,"error reading IO_INFILE")
       endif

!pw emep1.2beta

       fi_EMEP = -32.0
       an_EMEP = 237.7316364 ! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
       xp_EMEP =  43.0
       yp_EMEP = 121.0

       if(abs(fi_EMEP-fi)+abs(an_EMEP-an)+ &
          abs(xp_EMEP-xp)+abs(yp_EMEP-yp) > 0.0001)then

!Convert from EMEP coordinates to present coordinates if the map 
!is not the EMEP map

          write(*,*)'Converting fields read from ',fname,' to present map'
          write(*,*)fi_EMEP,an_EMEP,xp_EMEP,yp_EMEP
          write(*,*)fi,an,xp,yp

          imaxout = GIMAX 
          jmaxout = GJMAX 
          xpout = xp -ISMBEG +1
          ypout = yp -JSMBEG +1
          fiout = fi
          anout = an
          imaxin = IILARDOM
          jmaxin = JJLARDOM
       call ij2ij(in_field,imaxin,jmaxin,out_field,imaxout,jmaxout, &
                   fi_EMEP,an_EMEP,xp_EMEP,yp_EMEP, &
                   fiout,anout,xpout,ypout)

! Shift array (for compatibility) 
       do j = JSMBEG,JSMBEG+GJMAX-1
          j0 = j-JSMBEG+1
          do i = ISMBEG,ISMBEG+GIMAX-1
             i0 = i - ISMBEG+1
             in_field(i,j) = out_field(i0,j0)
          enddo
       enddo

       endif

    endif !me==0
            
    call global2local(in_field,local_field,MSG_READ7                 &
                       ,1,IILARDOM,JJLARDOM,1,ISMBEG,JSMBEG)

  end subroutine ReadField_r

!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
   subroutine ReadField_i(IO_INFILE,fname,local_field)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character*20, intent(in) :: fname         ! File name
  integer, intent(out)     :: local_field(MAXLIMAX,MAXLJMAX)

  !/-- We need an array for the whole model domain

  integer :: in_field(IILARDOM,JJLARDOM)! Field to be read
  real :: in_field_r(IILARDOM,JJLARDOM)
  real :: out_field(GIMAX ,GJMAX ) !pw


  integer :: imaxin,jmaxin,imaxout,jmaxout,i0,j0 !pw
  real :: fi_EMEP,an_EMEP,xp_EMEP,yp_EMEP,fiout,anout,xpout,ypout !pw

   if (me==0)then
      call open_file(IO_INFILE,"r",fname,needed=.true.)
      if ( ios /= 0 )then         
         write(6,*) 'error in opening IO_INFILE', fname
         call gc_abort(me,NPROC," error opening in_field")
      endif
    endif !me==0

!In the old code the field was assumed to be in the correct order
    if (me == 0) then

        READFIELD : do n = 1, BIG
           read(IO_INFILE,*,iostat=ios) i,j,in_field(i,j)
           if (  i < 1 .or. i > IILARDOM  .or. &
                 j < 1 .or. j > JJLARDOM  ) then  
              write(6,*)'error in i,j index in IO_INFILE=',fname
              ios = 88  !! Set an error code
              exit READFIELD
           endif
           if ( ios /= 0 ) exit READFIELD
        enddo READFIELD

        close(IO_INFILE)

        if ( ios /= 0 )then             
           write(6,*) 'error in reading IO_INFILE', fname
           call gc_abort(me,NPROC,"error reading IO_INFILE")
        endif !ios
!pw emep1.2beta
       fi_EMEP = -32.0
       an_EMEP = 237.7316364 ! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
       xp_EMEP =  43.0
       yp_EMEP = 121.0

       if(abs(fi_EMEP-fi)+abs(an_EMEP-an)+ &
          abs(xp_EMEP-xp)+abs(yp_EMEP-yp) > 0.0001)then


!Convert from EMEP coordinates to present coordinates if the map 
!is not the EMEP map

          write(*,*)'Converting fields read from ',fname,' to present map'
          in_field_r = real(in_field)
          
          imaxout = GIMAX
          jmaxout = GJMAX
          xpout = xp -ISMBEG +1
          ypout = yp -JSMBEG +1
          fiout = fi
          anout = an
          imaxin = IILARDOM
          jmaxin = JJLARDOM
       call ij2ij(in_field_r,imaxin,jmaxin,out_field,imaxout,jmaxout, &
                   fi_EMEP,an_EMEP,xp_EMEP,yp_EMEP, &
                   fiout,anout,xpout,ypout)

! Shift array (for compatibility)
       do j = JSMBEG,JSMBEG+GJMAX-1
          j0 = j-JSMBEG+1
          do i = ISMBEG,ISMBEG+GIMAX-1
             i0 = i - ISMBEG+1
             in_field(i,j) = nint(out_field(i0,j0))
          enddo
       enddo

       endif


    endif !me==0

    call global2local_int(in_field,local_field,MSG_READ5               &
               ,IILARDOM,JJLARDOM,1,ISMBEG,JSMBEG)
 
  end subroutine ReadField_i

end module ReadField_ml
