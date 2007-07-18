                         module Out_restri_ml
!    The coordinates for output on a restricted dommain are set in 
!    My_Outputs_ml. The restricted output may serve as boundary 
!    concentrations for computations with finer scale models.
!
!    To avoid this output:
!        ,ISPEC_OUTEND = -130  &     ! Set negative in My_Outputs_ml
!        ,JSPEC_OUTEND = -80         ! Set negative in My_Outputs_ml
!
!    If the coordinates of the restricted dommain are negative 
!    ( upper bound<than lower bound ) or if the restricted dommain is 
!    extending outside the computing dommain, then NO RESTRICTED OUTPUT 
!    IS CREATED!
!    
!
      use ModelConstants_ml, only: NPROC
      use Par_ml  , only : IRUNBEG,JRUNBEG, me
      use My_Outputs_ml, only: ISPEC_OUTBEG,JSPEC_OUTBEG  &
                              ,ISPEC_OUTEND ,JSPEC_OUTEND

implicit none
private

public :: set_outrestri


!	derived parameters - i,j postions in the actual grid (the same 
!       grid as the meteorological fields read in from infield )

  integer, public, parameter ::  &
         ISPEC_OUTBEG_ACT = ISPEC_OUTBEG - IRUNBEG + 1 &
        ,JSPEC_OUTBEG_ACT = JSPEC_OUTBEG - JRUNBEG + 1 &
        ,ISPEC_OUTEND_ACT = ISPEC_OUTEND - IRUNBEG + 1  &
        ,JSPEC_OUTEND_ACT = JSPEC_OUTEND - JRUNBEG + 1

!	output requested/possible, will be checked here

 logical, public, save :: to_out_restri

!	length in both directions+2D

  integer, public, parameter ::  &
         ILEN_RESTRI=ISPEC_OUTEND-ISPEC_OUTBEG+1     &		
        ,JLEN_RESTRI = JSPEC_OUTEND-JSPEC_OUTBEG+1

  integer, public, parameter ::  &
         MFSIZE_RESTRI=ILEN_RESTRI*JLEN_RESTRI

!	tg are the adresses in the output array

integer, public, save, dimension(0:NPROC-1) ::  &
      tlimax_restri, tgi0_restri    &
     ,tljmax_restri, tgj0_restri

!	tl are the adresses on the processors (>=1,<=tlimax(n))

integer, public, save, dimension(0:NPROC-1) ::  &
     tli0_restri, tli1_restri    &
    , tlj0_restri, tlj1_restri   &
    , tldim_restri

contains
subroutine set_outrestri()

        use Par_ml  , only : MAXLIMAX,MAXLJMAX &
                            ,GIMAX,GJMAX,me   &
                            ,tgi0,tgj0,tgi1,tgj1,tlimax,tljmax
        implicit none
        integer n
	
        to_out_restri = .true.

!      check if the restricted domain fits into the actual computing grid

       if( (ISPEC_OUTBEG_ACT.lt.1) .or. (JSPEC_OUTBEG_ACT.lt.1)  &
        .or. (ISPEC_OUTEND_ACT.gt.GIMAX)    &
        .or. (JSPEC_OUTEND_ACT.gt.GJMAX)) then
          print *,'out_restri will not be done'
          print *,'restricted domain not inside computational domain'
          print *,'computational domain in E/W is',IRUNBEG,IRUNBEG+GIMAX-1
          print *,'computational domain in S/N is',JRUNBEG,JRUNBEG+GJMAX-1
          print *,'out-restri domain in E/W is',ISPEC_OUTBEG,ISPEC_OUTEND
          print *,'out-restri domain in S/N is',JSPEC_OUTBEG,JSPEC_OUTEND
          to_out_restri = .false.
       endif

!      check, if upper bound>lower bound

       if((ISPEC_OUTBEG.gt.ISPEC_OUTEND) .or. &
          (JSPEC_OUTBEG.gt.JSPEC_OUTEND) )then
          if ( me == 0 ) then !u4
            print *,'out_restri will not be done'
            print *,'upper bound of out-restri domain less than lower bound'
            print *,'out-restri domain in E/W is',ISPEC_OUTBEG,ISPEC_OUTEND
            print *,'out-restri domain in S/N is',JSPEC_OUTBEG,JSPEC_OUTEND
          end if
          to_out_restri = .false.
        endif

!      assign

        do n=0,NPROC-1
          tgi0_restri(n) = tgi0(n) - ISPEC_OUTBEG_ACT + 1
          if(ISPEC_OUTBEG_ACT.gt.tgi0(n))tgi0_restri(n) = 1
          tgj0_restri(n) = tgj0(n) - JSPEC_OUTBEG_ACT + 1
          if(JSPEC_OUTBEG_ACT.gt.tgj0(n))tgj0_restri(n) = 1
          tli0_restri(n) = 1
          if(ISPEC_OUTBEG_ACT.gt.tgi0(n))tli0_restri(n) 	&
             = ISPEC_OUTBEG_ACT - tgi0(n) + 1
          tlj0_restri(n) = 1
          if(JSPEC_OUTBEG_ACT.gt.tgj0(n))tlj0_restri(n) 	&
             = JSPEC_OUTBEG_ACT - tgj0(n) + 1
          tli1_restri(n) = tgi1(n) - tgi0(n) + 1
          if(ISPEC_OUTEND_ACT.lt.tgi1(n))tli1_restri(n) 	&
             = ISPEC_OUTEND_ACT - tgi0(n) + 1
          tlj1_restri(n) = tgj1(n) - tgj0(n) + 1
          if(JSPEC_OUTEND_ACT.lt.tgj1(n))tlj1_restri(n) 	&
             = JSPEC_OUTEND_ACT - tgj0(n) + 1
          tlimax_restri(n) = tli1_restri(n) - tli0_restri(n) + 1
          tljmax_restri(n) = tlj1_restri(n) - tlj0_restri(n) + 1
          if((tlimax_restri(n).gt.0).and.(tlimax_restri(n).gt.0))  &
               tldim_restri(n) = tlimax_restri(n)*tljmax_restri(n)

!      test printout

        if(me == 0 .and. to_out_restri )then
          print *,'tgi',n,tgi0(n),tgi1(n),tlimax(n)
          print *,'tgi_restri',n,tgi0_restri(n),tlimax_restri(n)
          print *,'tli_restri',n,tli0_restri(n),tli1_restri(n)
          print *,'tgj',n,tgj0(n),tgj1(n),tljmax(n)
          print *,'tgj_restri',n,tgj0_restri(n),tljmax_restri(n)
          print *,'tlj_restri',n,tlj0_restri(n),tlj1_restri(n)
          print *,'tldim_restri',n,tldim_restri(n)
        endif
        enddo

end subroutine set_outrestri

end module Out_restri_ml
