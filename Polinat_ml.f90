module Polinat_ml

 use My_Outputs_ml,   only : NADV_FLIGHT, FLIGHT_ADV !ds - added
 use Chemfields_ml  , only : xn_adv
 use GenSpec_adv_ml
 use GridValues_ml , only : gl, gb
 use Io_ml,          only : IO_AIRCR
 use Met_ml,         only : z_bnd,z_mid
 use ModelConstants_ml , only : dt_advec,current_date,PPBINV,KMAX_BND
 use Par_ml   , only : gi0,gi1,gj0,gj1,ISMBEG,JSMBEG,me,NPROC
 implicit none
 private

 public polinat_init
 public polinat_in
 public polinat_out

 integer, private, save :: iimax, iii
 integer, private, save ::  fapos(2,999)
 real, private, save ::  kfalc(999), rhour(999)

 contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	subroutine polinat_init
	implicit none

	iii = 1
	rhour(1) = 1.
	rhour(2) = 0.

	end subroutine polinat_init

   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	subroutine polinat_in

        character*20 falc
	logical polexist

	integer ii,info

	if(current_date%month .eq. 6) then
	  if(current_date%day.eq.1 .or. current_date%day.eq.16 ) then
		if(me.eq.0)then
		  write(falc,fmt='(''tra9606'',i2.2,''.pos'')') 	&
			current_date%day
		  inquire(file=falc,exist=polexist)
		  write(6,*)'polinat exists?',polexist
		  if(.not.polexist)goto 912
		  open(IO_AIRCR,file=falc,status='unknown')
		  ii = 0
		  do while (.true.)
		    ii = ii + 1
		    read(IO_AIRCR,*,end=701) rhour(ii), fapos(1,ii), 	&
			fapos(2,ii), kfalc(ii)
		  end do
701		  continue
		  iimax = ii
		  rhour(iimax+1) = 0.
		  write(6,*) 'falcon positions  ',iimax
		  write(6,*) (rhour(ii),			&
			fapos(1,ii), fapos(2,ii),kfalc(ii),ii=1,5)
		  close(IO_AIRCR)
		  open(IO_AIRCR,file='aircraft.dat',position='append')
		  write(IO_AIRCR,*) 'month and day ',current_date%month	&
			,current_date%day
		  close(IO_AIRCR)
		endif
		iii = 1
!su	read on node 0

912		continue

!su	now distribute

		call gc_ibcast(760, 1, 0, NPROC, info, iimax)
		call gc_rbcast(761, iimax+1, 0, NPROC, info, rhour)
		call gc_rbcast(762, iimax, 0, NPROC, info, kfalc)
		call gc_ibcast(763, 2*iimax, 0, NPROC, info, fapos)

!su	all distributed
	  endif
	endif
!	if(current_date%month .eq. 7) then
!c
!		if(me.eq.0)then
!		  write(falc,fmt='(''tra9607'',i2.2,''.pos'')') 	&
!			current_date%day	      
!		  open(IO_AIRCR,file=falc,status='unknown')
!		  ii = 0
!		  do while (.true.)
!		    ii = ii + 1
!		    read(IO_AIRCR,*,end=702) rhour(ii), fapos(1,ii), 	&
!			fapos(2,ii), kfalc(ii)
!		  end do
!702		  continue
!		  iimax = ii
!		  rhour(iimax+1) = 0.
!		  write(6,*) 'falcon positions  ',iimax
!		  write(6,*) (rhour(ii),			&
!			fapos(1,ii), fapos(2,ii),kfalc(ii),ii=1,5)
!		  close(IO_AIRCR)
!		  open(IO_AIRCR,file='aircraft.dat',position='append')
!		  write(IO_AIRCR,*) 'month and day ',current_date%month	&
!			,current_date%day
!		  close(IO_AIRCR)
!		endif
!		iii = 1
!csu	read on node 0
!csu	now distribute

!		call gc_ibcast(760, 1, 0, NPROC, info, iimax)
!		call gc_rbcast(761, iimax+1, 0, NPROC, info, rhour)
!		call gc_rbcast(762, iimax, 0, NPROC, info, kfalc)
!		call gc_ibcast(763, 2*iimax, 0, NPROC, info, fapos)
!csu	all distributed
!	endif

	return
	end subroutine polinat_in

   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	subroutine polinat_out

	real ttt, dtmil

	integer ii,jj,k,jjj,info
        integer :: i

!    POLINAT flight positions
!	if(me.eq.0) write(6,*) 'for tidsjekk',current_date%hour	&
!			,dt_advec,(rhour(ii),ii=1,5)

	dtmil = dt_advec/60./60.
	ttt = current_date%hour+current_date%seconds/3600.
	if (rhour(2).gt.rhour(1) 			&
		.and. ttt+dtmil .gt.rhour(1)) then
	  do jjj = 1,10
	    if(me.eq.0) write(6,*) 'inne i tidsjekk',	&
			iii,jjj,ttt,rhour(iii),rhour(iii+1),dt_advec,dtmil

	    if (ttt .gt. rhour(iii) 		&
			.and. ttt .lt. rhour(iii+1)) then
!su	we have to synchronise the processors, since for next jjj(iii)
!su	the aircraft can be on another processor !!!!

	      call gc_gsync(NPROC,info)

	      if(me.eq.0) write(6,*) 'inne i tidsjekk2'	&
			,fapos(1,iii),fapos(1,iii), ttt
	      if(gi0+ISMBEG-1.le.fapos(1,iii) .and. 	&
			gi1+ISMBEG-1.ge.fapos(1,iii) .and.	&
			gj0+JSMBEG-1.le.fapos(2,iii) .and. 	&
			gj1+JSMBEG-1.ge.fapos(2,iii)) then
		write(6,*) 'inne i tidsjekk3',me,kfalc(iii)
		ii = fapos(1,iii) - gi0-ISMBEG+2
		jj = fapos(2,iii) - gj0-JSMBEG+2
            do k = 1,KMAX_BND-1
              if(z_bnd(ii,jj,k).gt.kfalc(iii) .and.       &
                  z_bnd(ii,jj,k+1).lt.kfalc(iii)) then
		    write(6,*) 'inne i tidsjekk4',me,kfalc(iii)
		    open(IO_AIRCR,file='aircraft.dat'	&
				,position='append')
		!ds uni.1: remove IXADV_O3 and replace by loop over FLIGHT_ADV
		    write(IO_AIRCR,*) ttt				&
			,( xn_adv( FLIGHT_ADV(i),ii,jj,k)*PPBINV, i=1, NADV_FLIGHT)  &
                  ,k,z_mid(ii,jj,k), gb(ii,jj),gl(ii,jj)
		    close(IO_AIRCR)
		  end if
		end do
	      end if
	      iii = iii + 1
	    end if
	    ttt = ttt + dtmil*0.1
	  end do
	end if

	return
	end subroutine polinat_out
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Polinat_ml


