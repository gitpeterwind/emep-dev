c
	subroutine outrestri_int2(itag,ifile,ident
     &			,dim1,ind1,arrout,scal)
c
c	now we scale all output-arrays accordingly (to fit into interger*2)
c	consequently, we must not check again on output!!!
c	mark the corresponding lines by 'CCHECK'
c
c
c	output of a 2D distributed array in integer*2 unformatted
c	proceeds in following steps:
c	1. input is arrout, dim1 is 1st dimension, ind1 the corresponding
c		index of the species to be written
c		scal is the scaling factor independent of the value
c		, depending on the array scale (i.e. to come to PPB)
c	2. compute the max of the array, define the additional scale iscal
c		, written to ident(20), to fit into integer*2
c	3. scale with the resulting scale: scal*10**(-iscal), result itmp
c	4. convert to integer*2, result itmp2
c	5. send to node 0, build up the output array buff
c	6. write on host node to file - putflti2
c
c	rCVTerr prints max/min values of scaled array with/without
c	the additional scaling, called of does not fit into integer*2
c
c	lines marked 'ccheckmax': print out the max of the array
c		for a given specification ident(6) - to be
c		assigned to narr_tocheck
c		if in addition one wants to know max/min and the coordinates
c		of these points call CVTerr - remove  cprintmax
c
	use par_ml  , only : MAXLIMAX,MAXLJMAX
     &			,NPROC,me
	use out_restri_ml , only: ILEN_RESTRI,MFSIZE_RESTRI
     &		,tlimax_restri, tgi0_restri
     &		,tljmax_restri, tgj0_restri
     &		, tli0_restri, tli1_restri
     &		, tlj0_restri, tlj1_restri
     &		,tldim_restri
chf	use My_Runmode_ml ,only: stop_test
c
	implicit none
c
	integer itag,ifile,ident(20),ierr
	integer dim1,ind1
	real arrout(dim1,MAXLIMAX,MAXLJMAX),scal
c
	integer MAXPK4
	integer NUMHOR
	parameter (NUMHOR=MAXLIMAX*MAXLJMAX)
#ifdef _CRAY
	integer NUMHORP,NUMHOR4
	parameter (NUMHORP=NUMHOR+3)
	parameter (NUMHOR4 = NUMHORP/4)
	parameter (MAXPK4=(MFSIZE_RESTRI+3)/4)
	integer buff(MAXPK4)
	integer itmp2(NUMHOR4)
	integer cvt64to16fac
#else
	parameter (MAXPK4=4*((MFSIZE_RESTRI+3)/4))
	integer*2 buff(MAXPK4)
	integer*2 itmp2(NUMHOR)
	integer ida1,iloct
#endif
c
	real rtmp(NUMHOR)
	integer i,j,ida,d
	integer ierror
	integer info,iscal
c
	real arrmax,scaladd,scaltest,locscaltest
c
ccheckmax	integer narr_tocheck 
c
c
c
c	define the maximum
c
	arrmax = 0.
	if(tldim_restri(me).gt.0)then
	  ida = 1
	  do j = tlj0_restri(me),tlj1_restri(me)
	    do i = tli0_restri(me),tli1_restri(me)
	      rtmp(ida) = arrout(ind1,i,j)
	      arrmax = max(arrmax,rtmp(ida))
	      ida = ida+1
	    enddo
	  enddo
	  arrmax = scal*arrmax
	endif
c	
	call gc_rmax(1, NPROC, info, arrmax)
c
	  iscal = -999
	  if(arrmax.gt.1.e-4)then
	    scaltest = arrmax/3.2e4
	    locscaltest = log10(scaltest)
	    if(locscaltest.le.0.)iscal = -int(-locscaltest)
	    if(locscaltest.gt.0.)iscal = int(locscaltest)+1
	  endif
c
c	compute the additional scale
c
c	if(arrmax.lt.3.2e1) then
c	  iscal = -3
c	else if(arrmax.lt.3.2e2) then
c	  iscal = -2
c	else if(arrmax.lt.3.2e3) then
c	  iscal = -1
c	else if(arrmax.lt.3.2e4) then
c	  iscal = 0
c	else if(arrmax.lt.3.e5) then
c	  iscal = 1
c	else 
c	  iscal = 2
c	end if
c
	ident(20) = iscal

	if(iscal.eq.-999) goto 1234	! output will contain only zero
					! write now and only ident
c
c	scale
c
	scaladd = scal*10.**(-iscal)
c
ccheckmax	narr_tocheck = 
ccheckmax	if(ident(6).eq.narr_tocheck)then
ccheckmax	  if(me.eq.0) write(6,*) 'max of array with number ',ident(6)
ccheckmax     &			,'is ',arrmax
cprintmax	  call CVTerr(rtmp,tldim_restri(me),scal,ident
cprintmax     &		,me,tlimax_restri(me),tgi0_restri(me), tgj0_restri(me))
ccheckmax	endif
c
c	convert to integer*2, send to host, build output array buff
c
#ifdef _CRAY
C
	if(tldim_restri(me).gt.0)then
	  ierror = cvt64to16fac(tldim_restri(me),rtmp,scaladd,itmp2)
CCHECK	  if(ierror.gt.0)call CVTerr(rtmp,tldim_restri(me),scal,ident
CCHECK     &		,me,tlimax_restri(me),tgi0_restri(me), tgj0_restri(me))
	  if(me.ne.0)call gc_isend(itag, (tldim_restri(me)+3)/4, 0, info,
     &              itmp2, itmp2)
	endif
c
	if(me.eq.0)then
	  buff(MAXPK4) = 0
c
	  if(tldim_restri(0).gt.0)then
	    call iprepput(ILEN_RESTRI,tgi0_restri(0)
     &		,tlimax_restri(0)
     &		,tgj0_restri(0),tljmax_restri(0)
     &		,(tldim_restri(0)+3)/4,MAXPK4,itmp2,buff)
	  endif
	  do d = 1, NPROC-1
	    if(tldim_restri(d).gt.0)then
	      call gc_irecv(itag,(tldim_restri(d)+3)/4, d, info, 
     &           itmp2, itmp2)
	      call iprepput(ILEN_RESTRI,tgi0_restri(d)
     &		,tlimax_restri(d)
     &		,tgj0_restri(d),tljmax_restri(d)
     &		,(tldim_restri(d)+3)/4,MAXPK4,itmp2,buff)
	    endif
	  enddo
	endif
#else
CCHECK	ierror = 0
	if(me.eq.0)then
	  do i = MFSIZE_RESTRI+1,MAXPK4
	    buff(i) = 0
	  enddo
	  if(tldim_restri(0).gt.0)then
	    ida = 0
	    ida1 = 1
	    do j = 1,tljmax_restri(0)
              do i = 1,tlimax_restri(0)
	        iloct = nint(scaladd*rtmp(ida1))
	        buff(ida+i) = iloct
CCHECK	        if(iloct.ne.buff(ida+i))ierror=1
		ida1 = ida1+1
	      enddo
	      ida = ida + ILEN_RESTRI
	    enddo
CCHECK	    if(ierror.gt.0)call CVTerr(rtmp,tldim_restri(0),scal,ident
CCHECK     &		,me,tlimax_restri(me),tgi0_restri(me), tgj0_restri(me))
	  endif
	  do d = 1, NPROC-1
	    if(tldim_restri(d).gt.0)then
	      call gc_brecv(itag,2*tldim_restri(d), d, info, 
     &           itmp2, itmp2)
              ida = (tgj0_restri(d)-1)*ILEN_RESTRI + tgi0_restri(d) - 1
	      ida1 = 1
	      do j = 1,tljmax_restri(d)
	        do i = 1,tlimax_restri(d)
	          buff(ida+i) = itmp2(ida1)
		  ida1 = ida1+1
	        enddo
	        ida = ida + ILEN_RESTRI
	      enddo
	    endif
	  enddo
	else
	  if(tldim_restri(me).gt.0)then
	  do j = 1,tldim_restri(me)
	      iloct = nint(scaladd*rtmp(j))
	      itmp2(j) = iloct
CCHECK	      if(iloct.ne.itmp2(j))ierror=1
	  enddo
CCHECK	  if(ierror.gt.0)call CVTerr(rtmp,tldim_restri(me),scal,ident
CCHECK     &		,me,tlimax_restri(me),tgi0_restri(me), tgj0_restri(me))
	  call gc_bsend(itag, 2*tldim_restri(me), 0, info,
     &              itmp2, itmp2)
	  endif
	endif
#endif
c
c	write to output
c
	ierr = 0

	if(me.eq.0)then
	  call putflti2(ifile,ident,MAXPK4,buff
     &			,MFSIZE_RESTRI,ierr)
	  if (ierr .ne. 0) then
	    write(6,*)' write error in outrestri.'
	    close(ifile)
            call gc_abort(me,NPROC,"rror in outrestri")
	  endif
	endif
c
chf u2	call stop_test(.true.,me,NPROC,ierr,'error in outrestri')

1234	continue
c
	return
	end
c
