	subroutine global2local(gloarr,locarr,msnr
     &			,dim0,dimi,dimj,diml,ibeg,jbeg)
c
c	distribute a 'real' array gloarr among the processors to get locarr
c	the array may have maximum 4 dimensions (i.e. snapemis)
c	, where the dimensions to be distributed, are dimi,dimj
c	the input array gloarr may be already restricted or not
c
	use PAR_ML , only : 
     &		MAXLIMAX	! Maximum number of local points in longitude
     &		,MAXLJMAX	! Maximum number of local points in latitude
     &		,NPROC		! Actual total number of processors
     &		,tgi0		! start points for all processors in longitude
     &		,tgj0		! start points for all processors in latitude
     &		,tlimax		! number of points for all processors in longitude
     &		,tljmax		! number of points for all processors in latitude
     &		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
c
	implicit none
c
c	input
	integer msnr		! message number
	integer dim0		! first dimension, possibly = 1 (= NSECTORS for snapemis)
     &		,dimi		! dimension in longitude (= GIMAX or IILARDOM)
     &		,dimj		! dimension in latitude  (= GJMAX or JJLARDOM)
     &		,diml		! 4th dimension, possibly = 1 (= NCMAX for snapemis)
     &		,ibeg		! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM
     &		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        real gloarr(dim0,dimi,dimj,diml)		! Global array
c
c	output
        real locarr(dim0,MAXLIMAX,MAXLJMAX,diml)	! Local array
c
c	local
	integer info,i,j,d,n0,nl
c
	if (me .ne. 0) then
c
c	receive from host
c
	  call gc_rrecv(msnr, dim0*MAXLIMAX*MAXLJMAX*diml, 0, info,
     &              locarr, locarr)	  

c
	else ! me = 0
c
c	first send to the others
c
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  do n0 = 1,dim0
		    locarr(n0,i,j,nl) = gloarr(n0,tgi0(d)+ibeg-2+i
     &				,tgj0(d)+jbeg-2+j,nl)
		  enddo
		enddo
	      enddo
	    enddo
	    call gc_rsend(msnr,dim0*MAXLIMAX*MAXLJMAX*diml, d, info, 
     &                 locarr, locarr)
	  enddo
c
c	now assign processor 0 itself
c
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		do n0 = 1,dim0
		  locarr(n0,i,j,nl) = gloarr(n0,i+ibeg-1,j+jbeg-1,nl)
		enddo
	      enddo
	    enddo
	  enddo
c
	endif	! me=?
c
	return
	end
c
c
	subroutine global2local_int(gloarr,locarr,msnr
     &			,dimi,dimj,diml,ibeg,jbeg)
c
c	distribute an 'integer' array gloarr among the processors to get locarr
c	the array may have maximum 3 dimensions (i.e. landcode)
c	, where the dimensions to be distributed, are dimi,dimj
c	the input array gloarr may be already restricted or not
c
	use PAR_ML , only : 
     &		MAXLIMAX	! Maximum number of local points in longitude
     &		,MAXLJMAX	! Maximum number of local points in latitude
     &		,NPROC		! Actual total number of processors
     &		,tgi0		! start points for all processors in longitude
     &		,tgj0		! start points for all processors in latitude
     &		,tlimax		! number of points for all processors in longitude
     &		,tljmax		! number of points for all processors in latitude
     &		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
c
	implicit none
c
c	input
	integer msnr		! message number
	integer dimi		! dimension in longitude (= GIMAX or IILARDOM)
     &		,dimj		! dimension in latitude  (= GJMAX or JJLARDOM)
     &		,diml		! 3rd dimension, possibly = 1 (= NCMAX for landcode)
     &		,ibeg		! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM
     &		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        integer gloarr(dimi,dimj,diml)		! Global array
c
c	output
        integer locarr(MAXLIMAX,MAXLJMAX,diml)	! Local array
c
c	local
	integer info,i,j,d,nl
c
	if (me .ne. 0) then
c
c	receive from host
c
	  call gc_irecv(msnr, MAXLIMAX*MAXLJMAX*diml, 0, info,
     &              locarr, locarr)	  
c
	else ! me = 0
c
c	first send to the others
c
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i
     &				,tgj0(d)+jbeg-2+j,nl)
		enddo
	      enddo
	    enddo
	    call gc_isend(msnr,MAXLIMAX*MAXLJMAX*diml, d, info, 
     &                 locarr, locarr)
	  enddo
c
c	now assign processor 0 itself
c
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
	      enddo
	    enddo
	  enddo
c
	endif	! me = ?
c
	return
	end
c
	subroutine global2local_short(gloarr,locarr,msnr
     &			,dimi,dimj,diml,ibeg,jbeg)
c
c	distribute a 'short integer' (integer *2) array gloarr among the processors to get locarr
c	the array may have maximum 3 dimensions (i.e. landcode)
c	, where the dimensions to be distributed, are dimi,dimj
c	the input array gloarr may be already restricted or not
c

	use PAR_ML , only : 
     &		MAXLIMAX	! Maximum number of local points in longitude
     &		,MAXLJMAX	! Maximum number of local points in latitude
     &		,NPROC		! Actual total number of processors
     &		,tgi0		! start points for all processors in longitude
     &		,tgj0		! start points for all processors in latitude
     &		,tlimax		! number of points for all processors in longitude
     &		,tljmax		! number of points for all processors in latitude
     &		,me		! Address of processor, numbering starts at 0 in south-west corner of ground level
c
	implicit none

      INCLUDE 'mpif.h'
      INTEGER STATUS(MPI_STATUS_SIZE)
c
c	input
	integer msnr		! message number
	integer dimi		! dimension in longitude (= GIMAX or IILARDOM)
     &		,dimj		! dimension in latitude  (= GJMAX or JJLARDOM)
     &		,diml		! 3rd dimension, possibly = 1 (= NCMAX for landcode)
     &		,ibeg		! start point of the array in longitude, = 1 for dimi = GIMAX or = ISMBEG for dimi = IILARDOM
     &		,jbeg		! start point of the array in latitude, = 1 for dimj = GJMAX or = JSMBEG for dimj = JJLARDOM
        integer*2 gloarr(dimi,dimj,diml)		! Global array
c
c	output
        integer*2 locarr(MAXLIMAX,MAXLJMAX,diml)	! Local array
c
c	local
	integer info,i,j,d,nl
c
	if (me .ne. 0) then
c
c	receive from host
c
c	  call gc_irecv(msnr, MAXLIMAX*MAXLJMAX*diml, 0, info,
c     &              locarr, locarr)	  
!pw use directly MPI routine: gc should be removed anyway!
       CALL MPI_RECV(locarr, 2*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE, 0,msnr,
     *     MPI_COMM_WORLD, STATUS, INFO)
c
	else ! me = 0
c
c	first send to the others
c
	  do d = 1, NPROC-1
	    do nl = 1,diml
	      do j = 1, tljmax(d)
		do i = 1, tlimax(d)
		  locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i
     &				,tgj0(d)+jbeg-2+j,nl)
		enddo
	      enddo
	    enddo
!	    call gc_isend(msnr,MAXLIMAX*MAXLJMAX*diml, d, info, 
!     &                 locarr, locarr)
!pw use directly MPI routine: gc should be removed anyway!
	    
      CALL MPI_SEND(locarr, MAXLIMAX*MAXLJMAX*diml*2, MPI_BYTE, d, msnr,
     *     MPI_COMM_WORLD, info)

	  enddo
c
c	now assign processor 0 itself
c
	  do nl = 1,diml
	    do j = 1, tljmax(0)
	      do i = 1, tlimax(0)
		locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
	      enddo
	    enddo
	  enddo
c
	endif	! me = ?
c
	return
	end

