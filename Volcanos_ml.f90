module Volcanos_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!
!                    module Volcanos_ml
!
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! Reads Volcanoes.dat with i,j,k for Volcanos
!
!u3 - added "private", ordered alphabetically
!u3   and corrected Makefile

 use My_Emis_ml  ,only : QRCVOL,molwt

 use EmisDef_ml,   only : NSECTORS,ISNAP_NAT
 use GridValues_ml, only  : sigma_bnd
 use Io_ml,         only : ios, open_file, IO_VOLC
 use ModelConstants_ml, only : KMAX_BND,KMAX_MID,PT
 use Met_ml,only : ps,roa
 use Par_ml,     only : ISMBEG, JSMBEG,me,NPROC,&
                        gi0,gj0,gi1,gj1 
 use PhysicalConstants_ml,  only :  GRAV,  AVOG

 implicit none
 private


 !/* subroutines:

  public :: VolcGet
  public :: Set_Volc      ! u3 added
  public :: Scale_Volc    ! u3 added


!u3 - some re-fomatting to get lines < 78 characters
  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO

   integer, public, parameter  :: NMAX_VOLC = 3 ! Max no. volcanos
   integer, public, save, &
           dimension(NMAX_VOLC):: height_volc  ! Height of volc
   integer, public, save       :: nvolc = 0    ! No volcanoes on emislist.sox
   integer, public, save, &
           dimension(NMAX_VOLC):: i_volc, j_volc
   integer, public, save       :: volc_no = 0    ! volc number
   real, public, save          :: emis_volc(NMAX_VOLC) = 0.0 ! emis. of volc


   real, private, save, dimension(NMAX_VOLC) :: &
            rcemis_volc0 !varies every hour
   real, public, save, dimension(NMAX_VOLC)  :: &
            rcemis_volc  !varies every time-step (as ps changes)

    logical, private, parameter :: DEBUG_VULC = .false.

contains 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine VolcGet(height_volc)
!***********************************************************************
!    DESCRIPTION:
! Reads file Volcanoes.dat with i,j,k for Volcanos
! Returns to EmisSet with height of volcanos
! Input file:
! Volcanoes.dat
!***********************************************************************

  integer, intent(out), dimension(NMAX_VOLC) :: height_volc
  integer  nvolc_read,height,i,j          ! Local variables
  character*13 fname

     
     !u3 write(fname,fmt='(''Volcanos.dat'')')  

     fname = "Volcanoes.dat"   ! u3 - simpler than write(fname !
     !u3 write(*,*)'I open ',fname  !u3 - open_file will report this.
     ios=0

     !u3 call open_file(IO_VOLC,"r",fname,needed=.true.)
     call open_file(IO_VOLC,"r",fname,needed=.true.,skip=1)
     if ( ios /= 0 )then
        WRITE(*,*) 'MPI_ABORT:'," VOLCGET: STOP"
        call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
     endif
     height_volc(:)=0.0
     nvolc_read=0

      READVOLC: do ! Read in volcanos 
               read(IO_VOLC,*,iostat=ios) i,j,height
               write(*,*)'found i,j,heigh',i,j,height
               if ( ios /= 0 ) exit READVOLC

               i = i -ISMBEG+1     ! for RESTRICTED domain
               j = j -JSMBEG+1     ! for RESTRICTED domain

               !/** Set the volcano number to be the same 
               !    as in gridSOx

               do volc_no=1,nvolc
                  if ((i_volc(volc_no)==i).and.     &
                             (j_volc(volc_no)==j))then
                     height_volc(volc_no)=height
                     nvolc_read=nvolc_read+1
                     write(*,*)'Found volcano with height k=',height
                  endif
               enddo
      enddo READVOLC

      write(6,*) nvolc_read,' volcanos on volcanos.dat&
             &match volcanos on emislist.sox'
      write(6,*)nvolc,' volcanos found in emislist.sox'
      if (nvolc_read < nvolc)then
         WRITE(*,*) 'MPI_ABORT:',"Volc missing in Volcanos.dat"
         call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
      endif   
      close(IO_VOLC)
      ios=0

    end subroutine VolcGet
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Set_Volc

    implicit none

    !**Local variables
    integer            :: k,i,j
    real               :: ehlpcom0(KMAX_MID)

    rcemis_volc0(:) = 0.
    ehlpcom0(:) = 0.

    !/** Set volcano
    do volc_no=1,nvolc
       k=height_volc(volc_no)
       i=i_volc(volc_no)
       j=j_volc(volc_no)
       if ((i >= gi0).and.(i<=gi1).and.(j>= gj0).and.&
          (j<= gj1))then !on the correct processor
          ehlpcom0(k)=GRAV* 0.001*AVOG/ &
                      (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))
          rcemis_volc0(volc_no)=emis_volc(volc_no)*ehlpcom0(k)/molwt(QRCVOL)  
       endif
    enddo ! volc_no
   end subroutine Set_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine Scale_Volc

   implicit none
   integer i,j,k,i_l,j_l
   real ehlpcom

  !/** Scale volc emissions to get emissions in molecules/cm3/s
  do volc_no=1,nvolc
     k=height_volc(volc_no)
     i=i_volc(volc_no)
     j=j_volc(volc_no)
     if ((i >= gi0).and.(i<=gi1).and.(j>= gj0).and.&
         (j<= gj1))then !on the correct processor
        i_l = i -gi0 +1 !local i
        j_l = j- gj0 +1 !local j
        ehlpcom=roa(i_l,j_l,KMAX_BND-k,1)/(ps(i_l,j_l,1)-PT)
        rcemis_volc(volc_no)=rcemis_volc0(volc_no)*ehlpcom
        if ( DEBUG_VULC ) then
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)
        end if ! DEBUG
     endif
  enddo ! volc_no
  end subroutine Scale_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module Volcanos_ml








