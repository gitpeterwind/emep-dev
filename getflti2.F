	subroutine getflti2(ifile,ident,itmp
     &		,ierr)
cfpp$ noconcur r
c
c  16 bit input and unpack
c
c  input:
c     mode:   0 = read field	is now hardcoded
c             1 = read field, skip fields until time > itime
c             2 = read field if field time = itime
c                 (otherwise next read starts at the same field)
c           100 = read field identification
c                 (next read starts at the same field)
c           101 = read field identification, skip fields
c                 until time > itime
c                 (next read starts at the same (last) field)
c           102 = read field identification, skip fields
c                 until time >= itime
c                 (next read starts at the same (last) field)
c           200 = scan rest of the file and read field with
c                 matching identification, specified identification
c                 input in ident(1:20) where -32767 means any value
c           201 = scan the whole file and read field with
c                 matching identification, specified identification
c                 input in ident(1:20) where -32767 means any value
c            -1 = clean up after a file is closed, and the same
c                 file unit no. is used for another file.
c     ifile: file unit no.
c     MFSIZEINPUT:  length of fdata (max field size)
c
c  output:
c     ident(20): field identification
c     fdata(..): field (unscaled, according to identification)
c     ierr = 0:  read o.k.
c            1:  read error
c            2:  read error, end_of_file
c
c
c  warning: using file unit no. (not file name) to identify
c           files when storing field identification.
c           if more than one file is opened with the same
c           unit, use 'call getflt(-1,...)' after closing
c           a file to avoid errors.
c
c  computer dependant i/o methodes for:
c              1) computer='cray'     (integer*2 not available)
c              2) computer='not.cray' (integer*2 used)
c
c  dnmi/fou  19.08.1993  anstein foss
c
c
c  modified by Peter Wind 11.03.2002:
c  uses nx=ident(10) (read from file) as first dimension of array, 
c  instead of IILARDOM. (only modified for not _CRAY)
c
	use Par_ml  , only : MAXLIMAX,MAXLJMAX,IILARDOM,JJLARDOM
     &			,MFSIZEINP,NPROC,ISMBEG,JSMBEG,tgi0,tgj0,me
     &			,MSG_INIT3
     &			,tlimax,tljmax,limax,ljmax
chf u2	use My_Runmode_ml , only : stop_test
	implicit none
c
	integer NUMHOR4,MAXPK4
	integer d,i,j,info
c  MAXPK4: max record length in cray 64 bit integer words
#ifdef _CRAY
	parameter (NUMHOR4=(MAXLIMAX*MAXLJMAX+3)/4)
	parameter (MAXPK4 = (IILARDOM*JJLARDOM+3)/4)
	integer   idpack(5)
	integer ipack(6+MAXPK4)
	integer itmp(6+NUMHOR4)
	integer cvt16to64fac,cvt16to64,ierror
#else
	parameter (NUMHOR4=MAXLIMAX*MAXLJMAX)
	parameter (MAXPK4 = IILARDOM*JJLARDOM)
	integer ida,itp
	integer*2 idpack(20)
	integer*2 ipack(21+MAXPK4)
	integer*2 itmp(21+NUMHOR4)
#endif
c
c  input/output
	integer   ifile,ident(20),ierr,iteserr
c
c
c  cray uses ipack as a standard length integer
c
c
c
	integer i,isave,nsave,ios,ierror,info
	integer nxin,nyin,nword,npack
c
	ierr = 0
	iteserr = 0
c
	if(me.eq.0)then
c
	  ipack(1) = 0
c
	  if(ifile.lt.1  .or. ifile.gt.99) then
	    write(6,*) ' **getflt** error ** ifile: ',ifile
	    ierr = 1
c	    call stop_all('GETFLTI2')
            call gc_abort(me,NPROC,"GETFLTI2")
	  endif
	endif

chf u2	call stop_test(.true.,me,NPROC,ierr,'GETFLTI2')
c
	if(me.eq.0)then
c        a) read field identification (one record)
c        b) read field data           (one record)
c
c
#ifdef _CRAY
	  read(ifile,iostat=ios,err=910,end=920) (idpack(i),i=1,5)
	  do i=1,5
	    ipack(i+1) = idpack(i)
	  enddo
c  convert (packed) integer*2 to cray integer
	  ierror = cvt16to64(20,idpack,ident)
#else
	  read(ifile,iostat=ios,err=910,end=920) (idpack(i),i=1,20)
	  do i=1,20
	    ident(i)=idpack(i)
	    ipack(i+1)=idpack(i)
	  end do
#endif
c
c  not using extra geometry identification (after field data)
	  if(ident(9).gt.999) ident(9)=ident(9)/1000
c
c
c
c
	  nxin=ident(10)
	  nyin=ident(11)
	  nword=nxin*nyin
c
	  if(nword.gt.MFSIZEINP) then
	    write(6,*) ' **getflt** field length too big',
     +                        ' (input buffer MFSIZEINP too small)'
	    write(6,*) ' **          MFSIZEINP = ',MFSIZEINP
	    write(6,*) ' ** ident: ',(ident(i),i=1,11)
	    write(6,*) ' **        ',(ident(i),i=12,20)
	    write(6,*) ' ** nx,ny,nx*ny: ',nxin,nyin,nword
	    ierr=1
	    iteserr = 1
	    goto 990
	  endif
c
#ifdef _CRAY
	  npack=(nword+3)/4
	  read(ifile,iostat=ios,err=910,end=920)(ipack(i),i=7,npack+6)
c  convert (packed) integer*2 to cray integer
#else
	  npack=nword
	  read(ifile,iostat=ios,err=910,end=920) (ipack(i),i=22,npack+21)
c
#endif
	  goto 990
c
c
910	  ierr=1
	  iteserr = 1
	  write(6,*) ' **getflt** read error. file,iostat: ',ifile,ios
	  goto 990
920	  ierr=2
	  goto 990
c
c
990	  continue
	endif
c	  if(ierr.eq.1)call stop_all('getflti2ex')
          if(ierr.eq.1)call gc_abort(me,NPROC,"getflti2ex")
chf u2	  call stop_test(.true.,me,NPROC,iteserr,'getflti2ex')

	if(me.eq.0)then

	if(ierr.eq.2)then
	    ipack(1) = -999
	    do d = 1, NPROC-1
#ifdef _CRAY
	      call gc_isend(MSG_INIT3,6+NUMHOR4,d, 
     &                    info, itmp,ipack)
#else
	      call gc_bsend(MSG_INIT3,2*(21+NUMHOR4),d, 
     &                    info, itmp,ipack)
#endif
	    enddo
	    close(ifile)
	    return
	endif
c
c
c..hj.. scaling is done within the getflti2 subroutine!!!!!!!!
c
	  if (ident(1) .ne. 88) then
	    write(6,*) 'ERROR IN INFIELD : produsent =',ident(1)
	  endif
c
	  if (ident(2) .ne. 1841.and. ident(2) .ne. 1600) then
	    write(6,*) 'ERROR IN INFIELD : grid =',ident(2)
	  endif
c
#ifdef _CRAY
 	  do i=1,6
	    itmp(i) = ipack(i)
	  enddo
	  do d = 1, NPROC-1
c
	    call iprepsend(IILARDOM,tgi0(d)+ISMBEG-1,MAXLIMAX
     &		,tgj0(d)+JSMBEG-1,MAXLJMAX
     &		,MAXPK4,NUMHOR4
     &		,ipack(7),itmp(7))
	    call gc_isend(MSG_INIT3,6+NUMHOR4,d, 
     &                    info, itmp,itmp)
	  enddo
	  call iprepsend(IILARDOM,ISMBEG,MAXLIMAX
     &		,JSMBEG,MAXLJMAX
     &		,MAXPK4,NUMHOR4
     &		,ipack(7),itmp(7))
#else
	  do i=1,21
	    itmp(i) = ipack(i)
	  enddo
	  do d = 1, NPROC-1
	    ida = 21
cpw	    itp = 21+(tgj0(d)+JSMBEG-2)*IILARDOM + tgi0(d)+ISMBEG-2
	    itp = 21+(tgj0(d)+JSMBEG-2)*ipack(11) + tgi0(d)+ISMBEG-2
	    do j = 1,tljmax(d)
	      do i = 1,tlimax(d)
	        itmp(ida+i) = ipack(itp+i)
	      enddo
	      ida = ida + MAXLIMAX
cpw	      itp = itp + IILARDOM
	      itp = itp + ipack(11)
	    enddo
	    call gc_bsend(MSG_INIT3,2*(21+NUMHOR4),d, 
     &                    info, itmp,itmp)
	  enddo
	  ida = 21
cpw	  itp = 21+(JSMBEG-1)*IILARDOM + ISMBEG-1
	  itp = 21+(JSMBEG-1)*ipack(11) + ISMBEG-1
	  do j = 1,ljmax
	    do i = 1,limax
	      itmp(ida+i) = ipack(itp+i)
	    enddo
	    ida = ida + MAXLIMAX
cpw	    itp = itp + IILARDOM
	    itp = itp + ipack(11)
	  enddo
#endif
c
c
	else		! me.ne.0, always receive to get itmp(1)
c
#ifdef _CRAY
c
	  call gc_irecv(MSG_INIT3,6+NUMHOR4,0, 
     &                 info, itmp, itmp)
#else
	  call gc_brecv(MSG_INIT3,2*(21+NUMHOR4),0, 
     &                 info, itmp, itmp)
#endif
	  if(itmp(1).eq.-999)then
	    ierr = 2
	    return
	  endif
#ifdef _CRAY
	  ierror = cvt16to64(20,itmp(2),ident)
c
#else
c
	  do i=1,20
	    ident(i) = itmp(i+1)
	  enddo
c
#endif
c
	endif
c
	return
	end
c
c
	subroutine getmetfield(ident,itmp,array)
c
	use Par_ml  , only : MAXLIMAX,MAXLJMAX
	implicit none
c
	integer NUMHOR4
	integer i
	integer ident
#ifdef _CRAY
	parameter (NUMHOR4=(MAXLIMAX*MAXLJMAX+3)/4)
	integer itmp(6+NUMHOR4)
	integer cvt16to64fac,ierror
#else
	parameter (NUMHOR4=MAXLIMAX*MAXLJMAX)
	integer*2 itmp(21+NUMHOR4)
#endif
	real scale
	real array(MAXLIMAX*MAXLJMAX)
c
	scale = 10.**ident
c
#ifdef _CRAY
c
	ierror = cvt16to64fac(MAXLIMAX*MAXLJMAX,itmp(7),array,scale)
c
#else
	do i = 1,MAXLIMAX*MAXLJMAX
	  array(i) = scale * itmp(i+21)
	enddo
#endif
c
c
	return
	end
c
c
#ifdef _CRAY
c
	INTEGER FUNCTION cvt16to64(NUM, ARR16, ARR64)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Convert 16 bit twos complement integers (integer*2) to 64 bit
C     *  twos complement integers (integer*8). The 16 bit integers are
C     *  assumed to be packed into 64 bit words.
C     *
C     * Input:
C     *  NUM       - number of items to be converted
C     *  ARR16     - array of packed 16 bit entities
C     *
C     * Output:
C     *  ARR64     - array of 64 bit sign extended entities
C     *  cvt16to64 - status code, always 0 for 16->64 bits expansion.
C     *
C     * Side effects:
C     *
C     ******************************************************************
cCDIR$ INTEGER=64
	IMPLICIT NONE
	INTEGER NUM, ARR16(1), ARR64(NUM)
	INTEGER I, J
csu	INTEGER MASK16, SIGN16, SIGNEXT64, IS
CSU      PARAMETER (MASK16    = X'000000000000FFFF',
CSU     $           SIGN16    = X'0000000000008000',
CSU     $           SIGNEXT64 = X'FFFFFFFFFFFF0000')
CSU      INTRINSIC AND, OR, SHIFTR
c
C---  unpack integer*2 into integer*8
c
	integer k,itmp,shfcnt,ibase
	intrinsic shiftl, shifta
c
	k = and(num,3)
	j = shiftr(num,2)
c
	if(num-k .ne. 0) then
	  ibase = 1
cdir$ unroll 8
	  do i =  1,j
	    itmp  = arr16(i)
	    ARR64(ibase+0) = shifta(itmp,48)
	    ARR64(ibase+1) = shifta(shiftl(itmp,16),48)
	    ARR64(ibase+2) = shifta(shiftl(itmp,32),48)
	    ARR64(ibase+3) = shifta(shiftl(itmp,48),48)
	    ibase = ibase + 4
	  enddo
	endif
c
	if(k .ne. 0) then
	  itmp  = arr16(j+1)
	  ARR64(ibase) = shifta(itmp,48)
	  shfcnt = 16
	  ibase = ibase + 1
	  do i=2,k
            ARR64(ibase) = shifta(shiftl(itmp,shfcnt),48)
            shfcnt = shfcnt + 16
            ibase = ibase + 1
	  enddo
	endif

	cvt16to64 = 0
	return
	END
c
c
	INTEGER FUNCTION cvt16to64fac(NUM, ARR16, ARR64,FAC)
C     ******************************************************************
C     * Purpose:
C     *
C     *  Convert 16 bit twos complement integers (integer*2) to 64 bit
C     *  twos complement integers (integer*8). The 16 bit integers are
C     *  assumed to be packed into 64 bit words.
C     *
C     * Input:
C     *  NUM       - number of items to be converted
C     *  ARR16     - array of packed 16 bit entities
C     *
C     * Output:
C     *  ARR64     - array of 64 bit sign extended entities
C     *  cvt16to64 - status code, always 0 for 16->64 bits expansion.
C     *
C     * Side effects:
C     *
C     ******************************************************************
cCDIR$ INTEGER=64
	IMPLICIT NONE
	INTEGER NUM, ARR16(1)
	REAL ARR64(NUM),fac
c
	INTEGER I, J
csu	INTEGER MASK16, SIGN16, SIGNEXT64, IS
CSU      PARAMETER (MASK16    = X'000000000000FFFF',
CSU     $           SIGN16    = X'0000000000008000',
CSU     $           SIGNEXT64 = X'FFFFFFFFFFFF0000')
CSU      INTRINSIC AND, OR, SHIFTR
c
C---  unpack integer*2 into integer*8
c
	integer k,itmp,shfcnt,ibase,jtmp
	intrinsic shiftl, shifta
c
	k = and(num,3)
	j = shiftr(num,2)
c
	if(num-k .ne. 0) then
	  ibase = 1
cdir$ unroll 8
	  do i =  1,j
	    itmp  = arr16(i)
	    jtmp = shifta(itmp,48)
	    ARR64(ibase+0) = jtmp*fac
	    jtmp = shifta(shiftl(itmp,16),48)
	    ARR64(ibase+1) = jtmp*fac
	    jtmp = shifta(shiftl(itmp,32),48)
	    ARR64(ibase+2) = jtmp*fac
	    jtmp = shifta(shiftl(itmp,48),48)
	    ARR64(ibase+3) = jtmp*fac
	    ibase = ibase + 4
	  enddo
	endif
c
	if(k .ne. 0) then
	  itmp  = arr16(j+1)
	  jtmp = shifta(itmp,48)
	  ARR64(ibase) = jtmp*fac
	  shfcnt = 16
	  ibase = ibase + 1
	  do i=2,k
	    jtmp = shifta(shiftl(itmp,shfcnt),48)
	    ARR64(ibase) = jtmp*fac
	    shfcnt = shfcnt + 16
	    ibase = ibase + 1
	  enddo
	endif
c
	cvt16to64fac = 0
	return
	end
c
c
	subroutine iprepsend(idim,ibeg,iend,jbeg,jend
     &		,lenin,lenout,arrin,arrout)
	implicit none
c
	intrinsic shiftr
	intrinsic shiftl
	integer lenin,lenout
	integer arrin(lenin),arrout(lenout)
	integer idim,ibeg,iend,jbeg,jend
c
	integer llout,lloutpos,lljbeg,lljbeg4,lljbegrest,lloutbeg
	integer j,i
	integer kkbeg,kkl
	integer intmp,outtmp
#ifdef _CRAY
!
#else
	external irmycopstrit
#endif
c
	lljbeg = ibeg+(jbeg-1)*idim+3
	lloutbeg = 4
c
	do j=1,jend
	  lljbeg4 = lljbeg/4
	  lljbegrest = lljbeg-lljbeg4*4
	  llout = lloutbeg/4
	  lloutpos = lloutbeg-llout*4
	  if(lloutpos.ne.lljbegrest)then
#ifdef _CRAY
	    call rmycopstrite(arrout(llout),arrin(lljbeg4)
     &		,2*iend,2*lloutpos,2*lljbegrest)
#else
	    call irmycopstrit(arrout(llout),arrin(lljbeg4)
     &		,2*iend,2*lloutpos,2*lljbegrest)
#endif
	  else
	    if(lloutpos.ne.0)then
	      kkbeg = 4-lloutpos
	      kkl = (iend+3-kkbeg)/4
	      intmp = shiftr(shiftl(arrin(lljbeg4),16*lloutpos)
     &			,16*lloutpos)
	      outtmp = shiftl(shiftr(arrout(llout),16*kkbeg)
     &			,16*kkbeg)
	      arrout(llout)=or(outtmp,intmp)
	      llout = llout+1
	      lljbeg4 = lljbeg4+1
	    else
	      kkl = (iend+3)/4
	    endif
	    do i=0,kkl-1
	      arrout(llout+i)= arrin(lljbeg4+i)
	    enddo
	  endif
	  lljbeg = lljbeg+idim
	  lloutbeg = lloutbeg+iend
	enddo
c
	return
	end
#endif
