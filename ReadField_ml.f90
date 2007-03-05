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
  character*50 :: errmsg 

  !/ rv2_1_6, ds 17/2/2005:
  !  Initialisation to zero added, so now we do not need to input an array which
  !  covers the whole domain. 
  !
  ! Also, print used instead of write for error messages, and
  ! errmsg more explicit now.

  real :: in_field(IILARDOM,JJLARDOM)! Field to be read

  in_field(:,:)    = 0.0       ! Initialise - ds, 15/1/2005      
  local_field(:,:) = 0.0       ! Initialise - ds, 15/1/2005      
  errmsg = "ok"

    if (me==0)then
       call open_file(IO_INFILE,"r",fname,needed=.true.)
          if ( ios /= 0 )then
          print *, 'error in opening IO_INFILE ', fname, ios
            WRITE(*,*) 'MPI_ABORT: ', "newmonth:error opening in_field" 
            call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
          endif
    endif !me==0


    if (me == 0) then

          READFIELD : do n = 1, BIG
             read(IO_INFILE,*,iostat=ios) i,j,in_field(i,j)
            if ( ios /= 0 ) exit READFIELD
             if (  i < 1 .or. i > IILARDOM  .or. &
                   j < 1 .or. j > JJLARDOM  ) then  
                  errmsg = "error in i,j index in IO_INFILE="!!! ,fname, i,j
                  ios = 88  !! Set an error code
                  exit READFIELD
             endif
          enddo READFIELD

       close(IO_INFILE)
       if ( errmsg /= "ok" ) then      
           print *, 'error in reading IO_INFILE', fname
           print *,  errmsg
             WRITE(*,*) 'MPI_ABORT: ', "errorreading IO_INFILE" 
             call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
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
  character*50 :: errmsg 

  !/-- We need an array for the whole model domain
  !/-- We *do not* need an array for the whole model domain!
  !    - changed ds, 15/1/2005

  integer :: in_field(IILARDOM,JJLARDOM)! Field to be read
  real :: in_field_r(IILARDOM,JJLARDOM)

  errmsg = "ok"

   if (me==0)then
      call open_file(IO_INFILE,"r",fname,needed=.true.)
      if ( ios /= 0 )then         
         print *, 'error in opening IO_INFILE', fname
           WRITE(*,*) 'MPI_ABORT: ', "error opening in_field" 
           call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
      endif
    endif !me==0

!In the old code the field was assumed to be in the correct order
    if (me == 0) then

        READFIELD : do n = 1, BIG
           read(IO_INFILE,*,iostat=ios) i,j,in_field(i,j)
           if ( ios /= 0 ) exit READFIELD
           if (  i < 1 .or. i > IILARDOM  .or. &
                 j < 1 .or. j > JJLARDOM  ) then  
              errmsg = "error in i,j index in IO_INFILE=" // fname
              ios = 88  !! Set an error code
              exit READFIELD
           endif
           if ( ios /= 0 ) exit READFIELD
        enddo READFIELD

        close(IO_INFILE)

       !ds if ( ios /= 0 )then             
       if ( errmsg /= "ok" ) then      
           print *, 'error in reading IO_INFILE', fname
             WRITE(*,*) 'MPI_ABORT: ', "errorreading IO_INFILE" 
             call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
        endif !ios

    endif !me==0

    call global2local_int(in_field,local_field,MSG_READ5               &
               ,IILARDOM,JJLARDOM,1,ISMBEG,JSMBEG)
 
  end subroutine ReadField_i

end module ReadField_ml
