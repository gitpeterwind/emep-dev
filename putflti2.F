c
	subroutine outarr_int2(itag,ifile,ident
     &			,dim1,ind1,arrout,scal)
c
c	now we scale all output-arrays accordingly (to fit into interger*2)
c	consequently, we must not check again on output!!!
c	mark the corresponding lines by 'CCHECK'
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
c	CVTerr prints max/min values of scaled array with/without
c	the additional scaling, called of does not fit into integer*2
c
c	lines marked 'ccheckmax': print out the max of the array
c		for a given specification ident(6) - to be
c		assigned to narr_tocheck
c		if in addition one wants to know max/min and the coordinates
c		of these points call CVTerr - remove  cprintmax
c
	use Par_ml  , only : MAXLIMAX,MAXLJMAX,MFSIZEOUT
     &			,GIMAX,GJMAX,limax,ljmax,NPROC,me
     &			,tgi0,tgj0,tlimax,tljmax
chf	use My_Runmode_ml ,only: stop_test
c
	implicit none
c
	integer itag,ifile,ident(20)
	integer ierr
	integer dim1,ind1
	real arrout(dim1,MAXLIMAX,MAXLJMAX),scal
c
	integer MAXPK4
	integer NUMHOR
	parameter (NUMHOR=MAXLIMAX*MAXLJMAX)
#ifdef _CRAY
	integer NUMHORP,NUMHOR4
	parameter (MAXPK4=(MFSIZEOUT+3)/4)
	parameter (NUMHORP=NUMHOR+3)
	parameter (NUMHOR4 = NUMHORP/4)
	integer buff(MAXPK4)
	integer itmp2(NUMHOR4)
	integer cvt64to16fac
#else
	parameter (MAXPK4=4*((MFSIZEOUT+3)/4))
	integer*2 buff(MAXPK4)
	integer*2 itmp2(NUMHOR)
	integer iloct,ida1
#endif
c
	real rtmp(NUMHOR)
	integer i,j,ida,d
	integer ierror
	integer info,iscal
c
	real arrmax,arrmin,scaladd,scaltest,locscaltest
c
ccheckmax	integer narr_tocheck 
c
c
c	define the maximum
c
	arrmax = -1.e40
	arrmin = 1.e40
	ida = 1
	do j = 1,ljmax
	  do i = 1,limax
	    rtmp(ida) = arrout(ind1,i,j)
	    arrmax = max(arrmax,rtmp(ida))
	    arrmin = min(arrmin,rtmp(ida))
	    ida = ida+1
	  enddo
	enddo
c
	arrmax = max(arrmax,-arrmin)
	arrmax = scal*arrmax
c
	
	call gc_rmax(1, NPROC, info, arrmax)
c
	  iscal = 0
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
c	scale
c
	scaladd = scal*10.**(-iscal)
c
	ident(20) = iscal
c
ccheckmax	narr_tocheck = 
ccheckmax	if(ident(6).eq.narr_tocheck)then
ccheckmax	  if(me.eq.0) write(6,*) 'max of array with number ',ident(6)
ccheckmax     &			,'is ',arrmax
cprintmax	  call CVTerr(rtmp,limax*ljmax,scal,ident
cprintmax     &			,me,tlimax(me),tgi0(me), tgj0(me))
ccheckmax	endif
c
c	convert to integer*2, send to host, build output array buff
c
#ifdef _CRAY
C
	ierror = cvt64to16fac(limax*ljmax,rtmp,scaladd,itmp2)
CCHECK	if(ierror.gt.0)call CVTerr(rtmp,limax*ljmax,scal,ident
CCHECK     &			,me,tlimax(me),tgi0(me), tgj0(me))
	if(me.eq.0)then
	  buff(MAXPK4) = 0
c
	  call iprepput(GIMAX,tgi0(0),tlimax(0)
     &		,tgj0(0),tljmax(0)
     &		,(limax*ljmax+3)/4,MAXPK4,itmp2,buff)
	  do d = 1, NPROC-1
	    call gc_irecv(itag,(tlimax(d)*tljmax(d)+3)/4, d, info, 
     &           itmp2, itmp2)
	    call iprepput(GIMAX,tgi0(d),tlimax(d)
     &		,tgj0(d),tljmax(d)
     &		,(tlimax(d)*tljmax(d)+3)/4,MAXPK4,itmp2,buff)
	  enddo
	else
	  call gc_isend(itag, (limax*ljmax+3)/4, 0, info,
     &              itmp2, itmp2)
	endif
#else
CCHECK	ierror = 0
	if(me.eq.0)then
	  do i = GIMAX*GJMAX+1,MAXPK4
	    buff(i) = 0
	  enddo
	  ida = 0
	  ida1 = 1
	  do j = 1,ljmax
	    do i = 1,limax
	      iloct = nint(scaladd*rtmp(ida1))
	      buff(ida+i) = iloct
CCHECK	      if(iloct.ne.buff(ida+i))error=1
	      ida1 = ida1+1
	    enddo
	    ida = ida + GIMAX
	  enddo
CCHECK	  if(ierror.gt.0)call CVTerr(rtmp,limax*ljmax,scal,ident
CCHECK     &			,me,tlimax(me),tgi0(me), tgj0(me))
	  do d = 1, NPROC-1
	    call gc_brecv(itag,2*tlimax(d)*tljmax(d), d, info, 
     &           itmp2, itmp2)
	    ida = (tgj0(d)-1)*GIMAX + tgi0(d)-1
	    ida1 = 1
	    do j = 1,tljmax(d)
	      do i = 1,tlimax(d)
	        buff(ida+i) = itmp2(ida1)
		ida1 = ida1+1
	      enddo
	      ida = ida + GIMAX
	    enddo
	  enddo
	else
	  do j = 1,limax*ljmax
	      iloct = nint(scaladd*rtmp(j))
	      itmp2(j) = iloct
CCHECK	      if(iloct.ne.itmp2(j))ierror=1
	  enddo
CCHECK	  if(ierror.gt.0)call CVTerr(rtmp,limax*ljmax,scal,ident
CCHECK     &			,me,tlimax(me),tgi0(me), tgj0(me))
	  call gc_bsend(itag, 2*limax*ljmax, 0, info,
     &              itmp2, itmp2)
	endif
#endif
c
c	write to output
c
	ierr = 0
	if(me.eq.0)then
	  call putflti2(ifile,ident,MAXPK4,buff,MFSIZEOUT,ierr)
	  if (ierr .ne. 0) then
	    write(6,*)' write error in outsig.'
	    close(ifile)
            call gc_abort(me,NPROC,"error in outchem")
	  endif
	endif
c
chf u2	call stop_test(.true.,me,NPROC,ierr,'error in outchem')
c
	return
	end
c
c
CCHECK	subroutine CVTerr(rtmp,tldim,scal,ident
CCHECK     &			,me,tlimaxme,tgi0me, tgj0me)
CCHECKc
CCHECK	implicit none
CCHECKc
CCHECK	integer ident(20),tldim
CCHECK	real rtmp(tldim),scal
CCHECK	integer me,tlimaxme,tgi0me, tgj0me
CCHECKc
CCHECK	real sc
CCHECK	real scaladd
CCHECK	real fmino,fmaxo
CCHECK	integer intfmin,intfmax
CCHECK	integer ifmin1d,ifmax1d
CCHECK	integer ifmin,ifmax
CCHECK	integer jfmin,jfmax,i,j,iscale
CCHECK	real fmin,fmax
CCHECKc
CCHECK	iscale = ident(20)
CCHECK	sc=10.**(-iscale)
CCHECK	scaladd = sc*scal
CCHECK	fmin=rtmp(1)
CCHECK	fmax=rtmp(1)
CCHECK	ifmin1d = 1
CCHECK	ifmax1d = 1
CCHECK	do j=2,tldim
CCHECK	    if(rtmp(j).lt.fmin)then
CCHECK	      fmin=rtmp(j)
CCHECK	      ifmin1d = j
CCHECK	    endif
CCHECK	    if(rtmp(j).gt.fmax)then
CCHECK	      fmax=rtmp(j)
CCHECK	      ifmax1d = j
CCHECK	    endif
CCHECK	end do
CCHECK	fmino = fmin*scal
CCHECK	fmaxo = fmax*scal
CCHECK	intfmin = nint(scaladd*fmin)
CCHECK	intfmax = nint(scaladd*fmax)
CCHECKc
CCHECK	j = ifmin1d/tlimaxme
CCHECK	i = ifmin1d - j*tlimaxme
CCHECK	ifmin = tgi0me+i-1
CCHECK	jfmin = tgj0me+j
CCHECKc
CCHECK	j = ifmax1d/tlimaxme
CCHECK	i = ifmax1d - j*tlimaxme
CCHECK	ifmax = tgi0me+i-1
CCHECK	jfmax = tgj0me+j
CCHECKc
CCHECK	write(6,*) ' **putflt** res64TO16 overflow. processor = ',me
CCHECK	write(6,*) ' **putflt** res64TO16 overflow.  field ident:'
CCHECK	write(6,*) ' ** ',(ident(i),i=1,11)
CCHECK	write(6,*) ' ** ',(ident(i),i=12,20)
CCHECK	write(6,*) ' ** field min,max:  ',fmino,fmaxo
CCHECK	write(6,*) ' ** scaled min,max: ',intfmin,intfmax,'  scaling: ',sc
CCHECK	write(6,*) ' ** min in point', ifmin,jfmin
CCHECK	write(6,*) ' ** max in point', ifmax,jfmax
CCHECKc
CCHECK	return
CCHECK	end
CCHECKc
c
	subroutine putflti2(ifile,ident,
     $                  ldata,ipack,output_size,ierr)
cfpp$ noconcur r
c
c  16 bit pack and output
c
c  input:
c     ifile: file unit
c     ident: field identification, ident(20)
c     lgeom: length of geometry section
c     igeom: geometry section, igeom(lgeom)
c     ldata: length of the field
c            (field size given in identification used when
c             writing the field)
c     fdata: the field, unscaled
c     scale: additional scaling of the field before output,
c            usually 1. (this is in addition to the scaling given
c            in the identification, for changing the basic unit)
c
c  output:
c     ierr = 0:  write o.k.
c            1:  write error
c
c  computer dependant i/o methodes for:
c              1) computer='cray'     (integer*2 not available)
c              2) computer='not.cray' (integer*2 used)
c
c  dnmi/fou  12.04.1992  anstein foss
c  dnmi/fou  19.08.1993  anstein foss
c
	implicit none
c
c  input/output
	integer   ifile,ident(20),ldata,ierr,output_size
c
c
c  cray uses ipack as a standard length integer
#ifdef _CRAY      
	integer cvt64to16
	integer   idpack(5),ipack(ldata)
#else
	integer*2 idpack(20),ipack(ldata)
#endif
	integer nxout,nyout,nword,i,nword4
	integer ierror
	integer lpack,ios,iloct,lidpack
c
c
	nxout=ident(10)
	nyout=ident(11)
	nword=nxout*nyout
c
	ierr=0
c
	if(nword.gt.output_size) then
	  write(6,*) ' **putflt** field length too big',
     +                        ' (buffer too small)'
	  write(6,*) ' **    output size = ',output_size
	  write(6,*) ' ** ident: ',(ident(i),i=1,11)
	  write(6,*) ' **        ',(ident(i),i=12,20)
	  write(6,*) ' ** nx,ny,nx*ny: ',nxout,nyout,nword
	  write(6,*) ' ** nword: ',nword
	  ierr=1
	  goto 990
	endif
c
c  data already scaled and 'packed'
c

c  always make output readable on cray (whole number of 64 bit words)
	nword4=((nword+3)/4)*4
c  zero padding
	ierror = 0
#ifdef _CRAY
c  convert cray integer to (packed) integer*2 ... identification
	ierror = cvt64to16(ident,idpack)
	lpack=nword4/4
	lidpack = 5
#else
	do i=1,20
	  idpack(i)=ident(i)
CCHECK	  iloct = int(idpack(i))
CCHECK	  if(iloct.ne.ident(i))ierror=1
	end do
	lpack=nword4
	lidpack = 20
#endif
CCHECK	if(ierror.ne.0) then
CCHECK	  write(6,*) ' **putflt** CVT64TO16 overflow.  ident:'
CCHECK	  write(6,*) ' ** ',(ident(i),i=1,11)
CCHECK	  write(6,*) ' ** ',(ident(i),i=12,20)
CCHECK	end if
c
c
c  write field identification and field data
	write(ifile,iostat=ios,err=90) (idpack(i),i=1,lidpack)
c	write now only, if iscale makes sense (cf. put_restri)
	if(ident(20).ne.-999)
     &		write(ifile,iostat=ios,err=90) (ipack(i),i=1,lpack)
c
	goto 990
c
90	ierr=1
	write(6,*) ' **putflt** write error. file,iostat: ',ifile,ios
c
990	continue
c
	return
	end
#ifdef _CRAY
	INTEGER FUNCTION cvt64to16fac(NUM, ARR64,scaladd, ARR16)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Convert 64 bit twos complement integers (integer*8) to 16 bit
C     *  twos complement integers (integer*2). The 16 bit integers are
C     *  assumed to be packed into 64 bit words.
C     *
C     * Input:
C     *  NUM       - number of items to be converted
C     *  ARR64     - array of 64 bit sign extended entities
C     *
C     * Output:
C     *  ARR16     - array of packed 16 bit entities
C     *  cvt64to16 - status code, 0 for successfull translation, count>0
C     *              if count items over(under)flowed.
C     *
C     * Side effects:
C     *
C     ******************************************************************
ccdir$ INTEGER=64
	IMPLICIT NONE
	INTEGER NUM, ARR16((NUM+3)/4)
	real ARR64(NUM),scaladd
c
************************************************************************
	INTEGER MININT16, MAXINT16, LIM16, I, J, ITMP
CCHECK	integer OFLOW
	PARAMETER (MININT16  = X'FFFFFFFFFFFF8000',
     $           MAXINT16  = X'0000000000007FFF')
	INTRINSIC OR, SHIFTL, SHIFTR
	intrinsic shifta
c
	integer k,jtmp,itmp1,itmp2,itmp3,shfcnt,ibase
	integer jtmp1,jtmp2,jtmp3,itmp4
c
C---  statement function to limit argument to within integer*2 range
	LIM16(I) = MIN(MAX(I,MININT16),MAXINT16)
c
CCHECK	oflow = 0
	k = and(num,3)
c
	if(num-k .ne. 0) then
	  ibase = 1
	  do i =  1,num-k,4
	    itmp   = nint(scaladd*arr64(i))
	    jtmp = shiftl(itmp,48)
CCHECK	    if(shifta(jtmp,48).ne.itmp) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp = shiftl(lim16(itmp),48)
CCHECK	    endif
	    itmp1  = nint(scaladd*arr64(i+1))
	    jtmp1 = shiftl(itmp1,48)
CCHECK	    if(shifta(jtmp1,48).ne.itmp1) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp1 = shiftl(lim16(itmp1),48)
CCHECK	    endif
	    itmp2  = nint(scaladd*arr64(i+2))
	    jtmp2 = shiftl(itmp2,48)
CCHECK	    if(shifta(jtmp2,48).ne.itmp2) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp2 = shiftl(lim16(itmp2),48)
CCHECK	    endif
	    itmp3  = nint(scaladd*arr64(i+3))
	    jtmp3 = shiftl(itmp3,48)
CCHECK	    if(shifta(jtmp3,48).ne.itmp3) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp3 = shiftl(lim16(itmp3),48)
CCHECK	    endif
	    itmp = or(jtmp,shiftr(jtmp1,16))
	    itmp4 = or(shiftr(jtmp2,32),shiftr(jtmp3,48))
csu            itmp = or(itmp,itmp4)
csu            arr16(ibase) = itmp
	    arr16(ibase) = or(itmp,itmp4)
	    ibase = ibase + 1
	  enddo
	endif

	if(k .ne. 0) then
	  itmp = nint(scaladd*arr64(num-k+1))
	  itmp1 = shiftl(itmp,48)
CCHECK	  if(shifta(itmp1,48).ne.itmp) then
CCHECK	    oflow = oflow + 1
CCHECK	    itmp1 = shiftl(lim16(itmp),48)
CCHECK	  endif
	  shfcnt = 16
	  do i=2,k
	    itmp = nint(scaladd*arr64(num-k+i))
	    jtmp = shiftl(itmp,48)
CCHECK	    if(shifta(jtmp,48).ne.itmp) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp = shiftl(lim16(itmp),48)
CCHECK	    endif
	    itmp1 = or(itmp1,shiftr(jtmp,shfcnt))
	    shfcnt = shfcnt + 16
	  enddo
	  arr16(ibase) = itmp1
	endif

	cvt64to16fac = 0
CCHECK	cvt64to16fac = OFLOW
	return
	END
c
c	now num is always=20
	INTEGER FUNCTION cvt64to16( ARR64, ARR16)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Convert 64 bit twos complement integers (integer*8) to 16 bit
C     *  twos complement integers (integer*2). The 16 bit integers are
C     *  assumed to be packed into 64 bit words.
C     *
C     * Input:
C     *  NUM       - number of items to be converted
C     *  ARR64     - array of 64 bit sign extended entities
C     *
C     * Output:
C     *  ARR16     - array of packed 16 bit entities
C     *  cvt64to16 - status code, 0 for successfull translation, count>0
C     *              if count items over(under)flowed.
C     *
C     * Side effects:
C     *
C     ******************************************************************
ccdir$ INTEGER=64
	IMPLICIT NONE
	INTEGER ARR16(5), ARR64(20)
c
************************************************************************
	INTEGER MININT16, MAXINT16, LIM16, I, J, ITMP
CCHECK	integer OFLOW
	PARAMETER (MININT16  = X'FFFFFFFFFFFF8000',
     $           MAXINT16  = X'0000000000007FFF')
	INTRINSIC OR, SHIFTL, SHIFTR
	intrinsic shifta
c
	integer k,jtmp,itmp1,itmp2,itmp3,shfcnt
	integer jtmp1,jtmp2,jtmp3,itmp4
c
C---  statement function to limit argument to within integer*2 range
	LIM16(I) = MIN(MAX(I,MININT16),MAXINT16)
c
CCHECK	oflow = 0
c
	  do i =  1,5
	    itmp   = arr64(4*(i-1)+1)
	    jtmp = shiftl(itmp,48)
CCHECK	    if(shifta(jtmp,48).ne.itmp) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp = shiftl(lim16(itmp),48)
CCHECK	    endif
	    itmp1  = arr64(4*(i-1)+2)
	    jtmp1 = shiftl(itmp1,48)
CCHECK	    if(shifta(jtmp1,48).ne.itmp1) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp1 = shiftl(lim16(itmp1),48)
CCHECK	    endif
	    itmp2  = arr64(4*(i-1)+3)
	    jtmp2 = shiftl(itmp2,48)
CCHECK	    if(shifta(jtmp2,48).ne.itmp2) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp2 = shiftl(lim16(itmp2),48)
CCHECK	    endif
	    itmp3  = arr64(4*(i-1)+4)
	    jtmp3 = shiftl(itmp3,48)
CCHECK	    if(shifta(jtmp3,48).ne.itmp3) then
CCHECK	      oflow = oflow + 1
CCHECK	      jtmp3 = shiftl(lim16(itmp3),48)
CCHECK	    endif
	    itmp = or(jtmp,shiftr(jtmp1,16))
	    itmp4 = or(shiftr(jtmp2,32),shiftr(jtmp3,48))
csu            itmp = or(itmp,itmp4)
csu            arr16(i) = itmp
	    arr16(i) = or(itmp,itmp4)
	  enddo


	cvt64to16 = 0
CCHECK	cvt64to16 = OFLOW
	return
	END
c
	subroutine iprepput(gidim,gibeg,gilen,gjbeg,gjlen
     &		,lenin,lenout,arrin,arrout)
	implicit none
c
	intrinsic shiftr
	intrinsic shiftl
	integer lenin,lenout
	integer arrin(lenin),arrout(lenout)
	integer gidim,gibeg,gilen,gjbeg,gjlen
c
	integer llin,llinpos,gllbeg,gllbeg4,gllbegrest,llinbeg
	integer j,i
	integer kkbeg,kkl,kkrest
	integer intmp,outtmp
#ifdef _CRAY
!
#else
	external irmycopstrit
#endif
c
	gllbeg = gibeg+(gjbeg-1)*gidim+3
	llinbeg = 4
	do j=1,gjlen
	  gllbeg4 = gllbeg/4
	  gllbegrest = gllbeg-gllbeg4*4
	  llin = llinbeg/4
	  llinpos = llinbeg-llin*4
	  if(llinpos.ne.gllbegrest)then
#ifdef _CRAY
	    call rmycopstrite(arrout(gllbeg4),arrin(llin)
     &		,2*gilen,2*gllbegrest,2*llinpos)
#else
	    call irmycopstrit(arrout(gllbeg4),arrin(llin)
     &		,2*gilen,2*gllbegrest,2*llinpos)
#endif
	  else
	    kkbeg = 0
	    if(llinpos.gt.0)then
	      kkbeg = 4-llinpos
	      intmp = shiftr(shiftl(arrin(llin),16*llinpos)
     &		,16*llinpos)
	      outtmp = shiftl(shiftr(arrout(gllbeg4),16*kkbeg)
     &		,16*kkbeg)
	      arrout(gllbeg4)=or(outtmp,intmp)
	      llin = llin+1
	      gllbeg4 = gllbeg4+1
	    endif
	    kkl = (gilen-kkbeg)/4
	    do i=0,kkl-1
	      arrout(gllbeg4+i)= arrin(llin+i)
	    enddo
	    kkrest = gilen-4*kkl-kkbeg
	    if(kkrest.gt.0)then
	      llin = llin+kkl
	      gllbeg4 = gllbeg4+kkl
	      kkbeg = 4-kkrest
	      intmp = shiftl(shiftr(arrin(llin),16*kkbeg),16*kkbeg)
	      outtmp = shiftr(shiftl(arrout(gllbeg4)
     &		,16*kkrest),16*kkrest)
	      arrout(gllbeg4)=or(outtmp,intmp)
	    endif
	  endif
	  gllbeg = gllbeg+gidim
	  llinbeg = llinbeg+gilen
	enddo
c
	return
	end
#endif
