module Volcanos_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!
!                    module Volcanos_ml
!
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !-----------------------------------------------------------------------!
  !  Preprocesses SOx emissions from volcanoes 
  !-----------------------------------------------------------------------!

 use CheckStop_ml,          only : CheckStop
 use My_Emis_ml,            only : QRCVOL,molwt
 use EmisDef_ml,            only : NSECTORS,ISNAP_NAT
 use GridValues_ml,         only : sigma_bnd, i_glob, j_glob
 use Io_ml,                 only : ios, open_file, check_file, IO_VOLC
 use ModelConstants_ml,     only : KMAX_BND,KMAX_MID,PT
 use Met_ml,                only : ps, roa
 use Par_ml,                only : ISMBEG, JSMBEG, me, NPROC, &
                                   li0,lj0,limax,ljmax
 use PhysicalConstants_ml,  only : GRAV, AVOG

 implicit none
 private


 !/* subroutines:

  public :: VolcGet
  public :: Set_Volc     
  public :: Scale_Volc   


  integer, public, parameter  :: NMAX_VOLC = 3  ! Max number of volcanoes
  integer, public, save       :: nvolc = 0    & ! No. grids with volcano 
                                                ! emissions in gridSOx
                                ,volc_no = 0 
  integer, public, save, dimension(NMAX_VOLC):: &
                        height_volc,          &  ! Height of volcanoes
                        i_volc, j_volc           ! Volcano's EMEP coordinates
  real, private, save, dimension(NMAX_VOLC)::   &
                        rcemis_volc0 ! Emissions part varying every hour
  real, public, save, dimension(NMAX_VOLC) ::   &
                        rcemis_volc, & ! Emissions part varying every time-step 
                        emis_volc = 0.0 ! Volcanoes' emissions

  logical, private, parameter :: DEBUG_VULC = .false.

contains 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine VolcGet(height_volc)
!***********************************************************************
 !-----------------------------------------------------------------------!
 ! Reads volcanoes' coorinates (i,j) and height level(height) 
 ! Returns to EmisSet with height of volcanos (height_volc)
 ! Input file: Volcanoes.dat
 !-----------------------------------------------------------------------!

  integer, intent(out), dimension(NMAX_VOLC) :: height_volc
  integer            :: nvolc_read,height,i,j          ! Local variables
  character (len=13) :: fname
  logical            :: fexist 
     
     fname = "Volcanoes.dat"  

     ios=0    ! Start with  assumed ok status

     call open_file(IO_VOLC,"r",fname,needed=.true.,skip=1)

     call CheckStop(ios,"VolcGet: problems with Volcanoes.dat ")


     height_volc(:)=0.0
     nvolc_read=0

     READVOLC: do
           read(IO_VOLC,*,iostat=ios) i,j,height

           write(*,*)'found i,j,heigh',i,j,height
           if ( ios /= 0 ) exit READVOLC

  !/** Read (i,j) are given for the full EMEP polar-stereographic domain
  !    Convert them to actual run domain
           i = i -ISMBEG+1    
           j = j -JSMBEG+1    

  !/** Set the volcano number to be the same as in emission data (gridSOx)

           do volc_no=1,nvolc
               if ((i_volc(volc_no)==i) .and. (j_volc(volc_no)==j)) then
                  height_volc(volc_no)=height
                  nvolc_read=nvolc_read+1
                  write(*,*)'Found volcano with height k=',height
               endif
           enddo
     enddo READVOLC

     write(6,*) nvolc_read,' volcanos on volcanos.dat &
             & match volcanos on emislist.sox'
     write(6,*) nvolc,' volcanos found in emislist.sox'

     call CheckStop(nvolc_read < nvolc, "Volc missing in Volcanos.dat")
  
     close(IO_VOLC)

    end subroutine VolcGet
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Set_Volc

 !-----------------------------------------------------------------------!
 ! Starts converting emission units from kg/m2/s to.... (hourly)
 !-----------------------------------------------------------------------!

    !**Local variables
    integer            :: k,i,j
    real               :: unit_conv1  

    rcemis_volc0(:) = 0.0
    unit_conv1      = 0.0

    !/** Set volcano
    do volc_no=1,nvolc
       k=height_volc(volc_no)
       i=i_volc(volc_no)
       j=j_volc(volc_no)
       if ( (i >= i_glob(li0)) .and. (i <= i_glob(limax) ) .and.  &
            (j >= j_glob(lj0)) .and. (j <= j_glob(ljmax)) ) then 

          unit_conv1 = GRAV* 0.001*AVOG/ &
                      (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))

          rcemis_volc0(volc_no) = emis_volc(volc_no)  &
                                * unit_conv1 / molwt(QRCVOL)
 
       endif
    enddo ! volc_no
   end subroutine Set_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine Scale_Volc

 !-----------------------------------------------------------------------!
 ! Finishing converting volcano emissions to molecules/cm3/s 
 ! (every advection timestep)
 !-----------------------------------------------------------------------!

   integer i,j,k,i_l,j_l
   real unit_conv2 


  do volc_no=1,nvolc

     k=height_volc(volc_no)
     i=i_volc(volc_no)
     j=j_volc(volc_no)

      if ( (i >= i_glob(li0)) .and. (i <= i_glob(limax) ) .and.  &
           (j >= j_glob(lj0)) .and. (j <= j_glob(ljmax)) ) then 

        i_l = i - i_glob(li0) + 1 !local i
        j_l = j - i_glob(lj0) + 1 !local j

        unit_conv2 = roa(i_l,j_l,KMAX_BND-k,1) / (ps(i_l,j_l,1)-PT)

        rcemis_volc(volc_no) = rcemis_volc0(volc_no) * unit_conv2

        if ( DEBUG_VULC ) &
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)

     endif
  enddo ! volc_no

  end subroutine Scale_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module Volcanos_ml
