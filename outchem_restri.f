c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     This file is the start of a parallel interface to the I/O routines
C     in Eulmod. Note that no parallel I/O is performed, it is only a
C     parallel interface to the I/O
C     
C     Content:
C     
C     outchem   -  output of air pollution fields in sigma levels
C     and accumulated wet and dry deposition at the surface
C
c	all the operations are now done in one call:
c		outarr_int2	cf. putflti2.F
c	here are only the calls and some presettings:
c		the specifications
c		the independent (characteristic for a given array)
c			scaling factor scale (to come i.e. to PPB)
c		total scale is then scal*10**(-iscal), iscal
c		is defined according to the maximum of the array
c		in outarr_int2 and written to ident(20), which is part
c		of the output fields
c		if the output array is some combination of other arrays,
c		then here this is assigned to rtmp	
C     
C----------------------------------------------------------------------------
C     $Id: outchem_restri.f,v 1.4 2002-09-13 08:55:39 mifads Exp $
C     Erik Berge, DNMI    Roar Skaalin, SINTEF Applied Mathematics
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
c     
	subroutine outchem_restri(outfilename)
c     
c     output of air pollution fields in sigma levels
c     and accumulated wet and dry deposition at the surface.
c     
c     the output parameters are: xn_adv,xn_shl,prodo3,accsu
c     
! uni - comment : atw is available also through GENOUT
c     
chf u2	use My_Runmode_ml ,only: stop_test
        use My_Outputs_ml ,only: 
     &                    ISPEC_OUTBEG, JSPEC_OUTBEG ! Start of restricted dom.
	use Par_ml, only:       me ,NPROC               ! processor number

	use Out_restri_ml , only : ILEN_RESTRI,JLEN_RESTRI
	use GenSpec_adv_ml, only:  NSPEC_ADV
c
	use ModelConstants_ml, only : PPBINV, identi
     &                  , KMAX_MID
      use GridValues_ml, only: sigma_mid, xp, yp
	use Io_ml,           only : IO_OUT          
	use Chemfields_ml  , only: xn_adv
	use GenSpec_adv_ml , only:  IXADV_O3, IXADV_NO, IXADV_NO2 
     >     ,IXADV_PAN, IXADV_HNO3, IXADV_MACR, IXADV_HCHO
     >     ,IXADV_CH3CHO, IXADV_C2H6, IXADV_NC4H10, IXADV_C2H4
     >     ,IXADV_C3H6, IXADV_OXYL, IXADV_ISOP, IXADV_H2O2
     >     ,IXADV_H2, IXADV_CO, IXADV_SO2, IXADV_SO4, IXADV_NITRATE
     >     ,IXADV_NH3, IXADV_AMSU, IXADV_GLYOX, IXADV_MGLYOX
     >     ,IXADV_MEK, IXADV_MVK, IXADV_MAL, IXADV_AMNI
	implicit none

c     
	character*30, intent(in) :: outfilename
c     
c----------------------------------------------------
c
c
c
c	local
	integer msgnr
	integer ident(20)
	integer ifile, i, k, n, ios, nn
	real scale
c     
c	open the output file
c
        ios = 0    
	if(me.eq.0)then
c     
	  ifile=IO_OUT
c     
	  open(ifile,file=outfilename,access='sequential',
     +        form='unformatted',status='unknown',iostat=ios)

	  if(ios.ne.0)then
	    write(6,*) 'can not open outfile ',ifile
	    write(6,*) 'filename : ',outfilename
            call gc_abort(me,NPROC,"error in outrestri")
	  endif

	endif

chf u2	call stop_test(.true.,me,NPROC,ios,'error in outrestri')
c     
	  ident(:) = identi(:)
c
c	do the output
c
	ident(4) = 6
	ident(5) = 2
        ident(8) = 0   ! jej suggestion
	ident(10) = ILEN_RESTRI
	ident(11) = JLEN_RESTRI

! Changes from jej, 10/5/01,  for position at the pole (xp, yp, read
! from infield).
! This is to make sure that the gridded data are located to the correct
! geographical location.

        ident(15) = 100* ( xp - (ISPEC_OUTBEG -1) )   ! position of N. pole
        ident(16) = 100* ( yp - (JSPEC_OUTBEG -1) )   ! ------" -----------
c     
	if(me.eq.0) write(6,*) 'outrestri i node 0',(ident(k),k=1,6)
c     
      do 1 k=1,KMAX_MID
	  ident(7)  = k
        ident(19) = sigma_mid(k)*10000.+0.5
c     
c     put the air-pollution variables
c     
	do n=1,NSPEC_ADV
	  ident(6)  = n+700
	  if(n .eq. IXADV_O3     .or.
     >       n .eq. IXADV_NO     .or. 
     >       n .eq. IXADV_NO2    .or. 
     >       n .eq. IXADV_PAN    .or. 
     >       n .eq. IXADV_HNO3   .or. 
     >       n .eq. IXADV_MACR   .or. 
     >       n .eq. IXADV_HCHO   .or. 
     >       n .eq. IXADV_CH3CHO .or. 
     >       n .eq. IXADV_C2H6   .or. 
     >       n .eq. IXADV_NC4H10 .or. 
     >       n .eq. IXADV_C2H4   .or.
     >       n .eq. IXADV_C3H6   .or.
     >       n .eq. IXADV_OXYL   .or.
     >       n .eq. IXADV_ISOP   .or.
     >       n .eq. IXADV_H2O2   .or.
     >       n .eq. IXADV_H2     .or.
     >       n .eq. IXADV_CO     .or.
     >       n .eq. IXADV_SO2    .or.
     >       n .eq. IXADV_SO4    .or.
     >       n .eq. IXADV_NITRATE.or.
     >       n .eq. IXADV_NH3    .or.
     >       n .eq. IXADV_AMSU   .or.
     >            n .eq. IXADV_AMNI     !  Ammonium nitrate
     >       .or. n .eq. IXADV_GLYOX    !  HCOHCO 
     >       .or. n .eq. IXADV_MGLYOX   !  H3COCHO = methyl glyoxal
     >       .or. n .eq. IXADV_MEK      !  CH3COC2H5 = methyl-ethyl-ketone
     >       .or. n .eq. IXADV_MVK      !  CH3C(=o)CH=CH2 = methyl-vinyl-ketone
     >       .or. n .eq. IXADV_MAL      !  CH3COCH=CHCHO = m..aldehyde 
     >       ) then
	  scale = PPBINV
	  msgnr = 100*n + k
	  call outrestri_int2(msgnr,ifile,ident
     &            ,NSPEC_ADV,n,xn_adv(1,1,1,k),scale)
c     
	  end if  ! species list
	enddo
1	continue
c
c     
	if (me.eq.0) close(ifile)
c     
c     
	return
c     
c     
	end
c
