	subroutine local2global(locarr,gloarr,msnr)
c
c	gather a 2D 'real' array from the processors at the host me=0
c
	use Par_ml , only : 
     &		MAXLIMAX	! Maximum number of local points in longitude
     &		,MAXLJMAX	! Maximum number of local points in latitude
     &		,GIMAX		! Number of global points in longitude
     &		,GJMAX		! Number of global points in latitude
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
c
	integer msnr			! message number
        real locarr(MAXLIMAX,MAXLJMAX)	! Local array
c
c	output
c
        real gloarr(GIMAX,GJMAX)	! Global array
c
c	local
	integer info,i,j,d
c
	if (me .ne. 0) then
c
c	send to host
c
	  call gc_rsend(msnr, MAXLIMAX*MAXLJMAX, 0, info,
     &              locarr, locarr)	  
c
	else ! me = 0
c
c	copy first local array
c
	  do j = 1, tljmax(0)
	    do i = 1, tlimax(0)
	      gloarr(i,j) = locarr(i,j)
	    enddo
	  enddo
c
c	now get from the others
c
	  do d = 1, NPROC-1
	    call gc_rrecv(msnr,MAXLIMAX*MAXLJMAX, d, info, 
     &                 locarr, locarr)
	    do j = 1, tljmax(d)
	      do i = 1, tlimax(d)
		gloarr(tgi0(d)-1+i,tgj0(d)-1+j)=locarr(i,j)
	      enddo
	    enddo
	  enddo
c
	endif	! me = ?
c
	return
	end
