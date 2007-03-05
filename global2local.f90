	subroutine global2local(gloarr,locarr,msnr&
     			,dim0,dimi,dimj,diml,ibeg,jbeg)
!
!	distribute a 'real' array gloarr among the processors to get locarr
!	the array may have maximum 4 dimensions (i.e. snapemis)
!	, where the dimensions to be distributed, are dimi,dimj
!	the input array gloarr may be already restricted or not
!
	use PAR_ML , only : &
     		MAXLIMAX&	! Maximum number of local points in longitude&
     		,MAXLJMAX&	! Maximum number of local points in latitude&
     		,NPROC	&	! Actual total number of processors&
     		,tgi0	&	! start points for all processors in longitude&
     		,tgj0	&	! start points for all processors in latitude&
     		,tlimax	&	! number of points for all processors in longitude&
     		,tljmax	&	! number of points for all processors in latitude&
     		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
!
	implicit none
!
        INCLUDE 'mpif.h'
        INTEGER STATUS(MPI_STATUS_SIZE),INFO

!	input
	integer msnr		! message number
	integer dim0&		! first dimension, possibly = 1 (= NSECTORS for snapemis)&
     		,dimi&		! dimension in longitude (= GIMAX or IILARDOM)&
     		,dimj&		! dimension in latitude  (= GJMAX or JJLARDOM)&
     		,diml&		! 4th dimension, possibly = 1 (= NCMAX for snapemis)&
     		,ibeg&		! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM&
     		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        real gloarr(dim0,dimi,dimj,diml)		! Global array
!
!	output
        real locarr(dim0,MAXLIMAX,MAXLJMAX,diml)	! Local array
!
!	local
	integer i,j,d,n0,nl
!
	if (me .ne. 0) then
!
!	receive from host
!
           CALL MPI_RECV( locarr, 8*dim0*MAXLIMAX*MAXLJMAX*diml, &
                MPI_BYTE,  0, msnr, MPI_COMM_WORLD, STATUS, INFO) 

!
	else ! me = 0
!
!	first send to the others
!
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  do n0 = 1,dim0
		    locarr(n0,i,j,nl) = gloarr(n0,tgi0(d)+ibeg-2+i&
     				,tgj0(d)+jbeg-2+j,nl)
		  enddo
		enddo
	      enddo
	    enddo
            CALL MPI_SEND(locarr,8*dim0*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE, &
                 d, msnr, MPI_COMM_WORLD, INFO) 
	  enddo
!
!	now assign processor 0 itself
!
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		do n0 = 1,dim0
		  locarr(n0,i,j,nl) = gloarr(n0,i+ibeg-1,j+jbeg-1,nl)
		enddo
	      enddo
	    enddo
	  enddo
!
	endif	! me=?
!
	return
	end
!
!
	subroutine global2local_int(gloarr,locarr,msnr&
     			,dimi,dimj,diml,ibeg,jbeg)
!
!	distribute an 'integer' array gloarr among the processors to get locarr
!	the array may have maximum 3 dimensions (i.e. landcode)
!	, where the dimensions to be distributed, are dimi,dimj
!	the input array gloarr may be already restricted or not
!
	use PAR_ML , only : &
     		MAXLIMAX&	! Maximum number of local points in longitude&
     		,MAXLJMAX&	! Maximum number of local points in latitude&
     		,NPROC	&	! Actual total number of processors&
     		,tgi0	&	! start points for all processors in longitude&
     		,tgj0	&	! start points for all processors in latitude&
     		,tlimax	&	! number of points for all processors in longitude&
     		,tljmax	&	! number of points for all processors in latitude&
     		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
!
	implicit none

        INCLUDE 'mpif.h'
        INTEGER STATUS(MPI_STATUS_SIZE),INFO
!
!	input
	integer msnr		! message number
	integer dimi	&	! dimension in longitude (= GIMAX or IILARDOM)&
     		,dimj&		! dimension in latitude  (= GJMAX or JJLARDOM)&
     		,diml	&	! 3rd dimension, possibly = 1 (= NCMAX for landcode)&
     		,ibeg	&	! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM&
     		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        integer gloarr(dimi,dimj,diml)		! Global array
!
!	output
        integer locarr(MAXLIMAX,MAXLJMAX,diml)	! Local array
!
!	local
	integer i,j,d,nl
!
	if (me .ne. 0) then
!
!	receive from host
!
	    CALL MPI_RECV(locarr, 4*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE,  0, &
	    msnr, MPI_COMM_WORLD,STATUS, INFO) 
!
	else ! me = 0
!
!	first send to the others
!
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i&
     				,tgj0(d)+jbeg-2+j,nl)
		enddo
	      enddo
	    enddo
	      CALL MPI_SEND( locarr, 4*MAXLIMAX*MAXLJMAX*diml, &
              MPI_BYTE, d, msnr, MPI_COMM_WORLD, INFO) 
	  enddo
!
!	now assign processor 0 itself
!
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
	      enddo
	    enddo
	  enddo
!
	endif	! me = ?
!
	return
	end
!
	subroutine global2local_short(gloarr,locarr,msnr&
     			,dimi,dimj,diml,ibeg,jbeg)
!
!	distribute a 'short integer' (integer *2) array gloarr among the processors to get locarr
!	the array may have maximum 3 dimensions (i.e. landcode)
!	, where the dimensions to be distributed, are dimi,dimj
!	the input array gloarr may be already restricted or not
!

	use PAR_ML , only : &
     		MAXLIMAX&	! Maximum number of local points in longitude&
     		,MAXLJMAX&	! Maximum number of local points in latitude&
     		,NPROC	&	! Actual total number of processors&
     		,tgi0	&	! start points for all processors in longitude&
     		,tgj0	&	! start points for all processors in latitude&
     		,tlimax	&	! number of points for all processors in longitude&
     		,tljmax	&	! number of points for all processors in latitude&
     		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
!
	implicit none

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE),INFO
!
!	input
	integer msnr		! message number
	integer dimi	&	! dimension in longitude (= GIMAX or IILARDOM)&
     		,dimj	&	! dimension in latitude  (= GJMAX or JJLARDOM)&
     		,diml	&	! 3rd dimension, possibly = 1 (= NCMAX for landcode)&
     		,ibeg	&	! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM&
     		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        integer*2 gloarr(dimi,dimj,diml)		! Global array
!
!	output
        integer*2 locarr(MAXLIMAX,MAXLJMAX,diml)	! Local array
!
!	local
	integer i,j,d,nl
!
	if (me .ne. 0) then
!
!	receive from host
!

       CALL MPI_RECV(locarr, 2*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE, 0,msnr,&
          MPI_COMM_WORLD, STATUS, INFO)
!
	else ! me = 0
!
!	first send to the others
!
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i&
     				,tgj0(d)+jbeg-2+j,nl)
		enddo
	      enddo
	    enddo

            CALL MPI_SEND(locarr, MAXLIMAX*MAXLJMAX*diml*2, MPI_BYTE, &
                 d, msnr,MPI_COMM_WORLD, info)

	  enddo
!
!	now assign processor 0 itself
!
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
	      enddo
	    enddo
	  enddo
!
	endif	! me = ?
!
	return
	end

