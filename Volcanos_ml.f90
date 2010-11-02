! <Volcanos_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 
module Volcanos_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!
!                    module Volcanos_ml
!
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !-----------------------------------------------------------------------!
  ! Processes SOx emission heights from volcanoes 
  ! In volcanos_ml, only the height of the volcanos are read from volcanos.dat,
  ! the emissions themselves comes from gridSOx, so this is to ensure that if
  ! suddenly e.g. Iceland starts to report volcano emissions and these are in
  ! gridSOx, the program will discover this.
  ! Note that nvolc is set according to the emission input from the gridSOx 
  ! files.  It counts the number of grids with the "country code" for volcanoes
  ! that are in the gridSOx file.
  !
  ! Note - EmisGet and other routines use "restricted" coords, which introduces
  ! some complications here. Hopefully we can tidy up one day.
  !-----------------------------------------------------------------------!

 use CheckStop_ml,          only : CheckStop
 use EmisDef_ml,            only : NSECTORS,ISNAP_NAT,VOLCANOES_LL
 use GridValues_ml,         only : GRIDWIDTH_M, xm2, sigma_bnd, i_fdom, j_fdom,i_local, j_local,lb2ij
 use Io_ml,                 only : ios, NO_FILE, open_file, check_file, IO_VOLC, Read_Headers,read_line
 use ModelConstants_ml,     only : KMAX_BND,KMAX_MID,PT, NPROC, MasterProc
 use MetFields_ml,          only : roa, ps
 use Par_ml,                only : IRUNBEG, JRUNBEG, me, li0,lj0,li1,lj1  &
                                  ,gi0, gi1, gj0, gj1 !TEST
 use PhysicalConstants_ml,  only : GRAV, AVOG
 use TimeDate_ml,           only : nydays  ! No. days per year
 use KeyValue_ml,           only : KeyVal, KeyValue, LENKEYVAL

 implicit none
 private


 !/* subroutines:

  public :: VolcGet
  public :: Set_Volc     
  public :: Scale_Volc   


  integer, public, parameter  :: NMAX_VOLC = 12  ! Max number of volcanoes
  integer, public, save       :: nvolc = 0    & ! No. grids with volcano 
                                                ! emissions in gridSOx
                                ,volc_no = 0 
  logical, public, save :: Volcanoes_found=.false.! true if Volcanoes found on this processor/subdomain

  integer, public, save, dimension(NMAX_VOLC):: &
                        height_volc,          &  ! Height of volcanoes
                        i_volc, j_volc           ! Volcano's EMEP coordinates
  real, private, save, dimension(NMAX_VOLC)::   &
                        rcemis_volc0 ! Emissions part varying every hour
  real, public, save, dimension(NMAX_VOLC) ::   &
                        rcemis_volc, & ! Emissions part varying every time-step 
                        emis_volc = 0.0 ! Volcanoes' emissions

  logical, private, parameter :: DEBUG_VULC = .true.

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO

contains 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine VolcGet(height_volc)
    !***********************************************************************
    !-----------------------------------------------------------------------!
    ! Reads volcanoes' coordinates (i,j) and height level(height) 
    ! Returns to EmisSet with height of volcanos (height_volc)
    ! Input file: Volcanoes.dat
    !-----------------------------------------------------------------------!

    integer, intent(out), dimension(NMAX_VOLC) :: height_volc
    integer            :: nvolc_read,height,i,j          ! Local variables
    character (len=23) :: fname
    logical            :: fexist 
    
    integer            ::nin,i_l,j_l
    real               ::lon,lat,value,xr,yr,conv,tonne_to_kgm2s

    character(len=70)               :: errmsg ! Message text
    character(len=20), dimension(4) :: Headers
    type(KeyVal), dimension(20)     :: KeyValues ! Info on units, coords, etc.
    integer                         :: NHeaders, NKeys
    character(len=80)               :: txtinput,s  ! Big enough to contain
                                                 ! one full input record
    ios=0    ! Start with  assumed ok status

    height_volc(:)=KMAX_MID
    nvolc_read=0
 
    if(.not.VOLCANOES_LL)then
       !Old style: data read from gridSOx and Volcanoes.dat

   fname = "Volcanoes.dat"  

    if (MasterProc)then
       call open_file(IO_VOLC,"r",fname,needed=.true.,skip=1)

       if(ios/=NO_FILE)then!"Volcanoes.dat" does exist
          call CheckStop(ios,"VolcGet: problems with Volcanoes.dat ")


          READVOLC: do
             read(IO_VOLC,*,iostat=ios) i,j,height

             if (DEBUG_VULC) write(*,"(a,2i5,f12.3,i9)") &
                   'VOLCIJ:found i,j,height,ios',i,j,height,ios
             if ( ios /= 0 ) exit READVOLC

             !/** Read (i,j) are given for the full EMEP polar-stereographic domain
             !    Convert them to actual run domain
             i = i -IRUNBEG+1    
             j = j -JRUNBEG+1    

             !/** Set the volcano number to be the same as in emission data (gridSOx)

             do volc_no=1,nvolc
                if ((i_volc(volc_no)==i) .and. (j_volc(volc_no)==j)) then
                   height_volc(volc_no)=height
                   nvolc_read=nvolc_read+1
                   if (DEBUG_VULC) write(*,*)'VOLCIJ:Found volcano with height k=',height
                endif
             enddo
          enddo READVOLC

          write(6,*) nvolc_read,' volcanos on volcanos.dat &
               & match volcanos on emislist.sox'
          write(6,*) nvolc,' volcanos found in emislist.sox'

          call CheckStop(nvolc_read < nvolc, "Volc missing in Volcanos.dat")


       endif
       close(IO_VOLC)
    endif

    CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,INFO)


    else ! use "VolcanoesLL.dat" and disregard volcanoes from gridSOx
       nvolc=0
       !read volcanoes from lon lat format file:
       fname = "VolcanoesLL.dat"  
       if (MasterProc) call open_file(IO_VOLC,"r",fname,needed=.true.)
       CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,INFO)

       if(ios/=0)then
          if (MasterProc) write(6,*) ' No volcanoes found in VolcanoesLL.dat'
          return
       endif

       call Read_Headers(IO_VOLC,errmsg,NHeaders,NKeys,Headers,Keyvalues)
       VOLCLOOP: do nin = 1, NMAX_VOLC+1
        call read_line(IO_VOLC,txtinput,ios)
        if ( ios /= 0 ) exit  ! End of file
        nvolc=nvolc+1
        call CheckStop(nvolc>NMAX_VOLC, "more Volcanoes in VolcanoesLL.dat than NMAX_VOLC")

        read(unit=txtinput,fmt=*) s,  lon,  lat, height_volc(nvolc), emis_volc(nvolc)
        call lb2ij(lon,lat,xr,yr)
        i_volc(nvolc)=nint(xr)
        j_volc(nvolc)=nint(yr)
        if (DEBUG_VULC.and.MasterProc) write(*,*)'VOLC:Found ',trim(s), &
             lon,  lat, i_volc(nvolc),j_volc(nvolc)&
             ,i_volc(nvolc)+IRUNBEG-1,j_volc(nvolc)+JRUNBEG-1, & 
              height_volc(nvolc), emis_volc(nvolc)

       end do VOLCLOOP
       if (MasterProc)close(IO_VOLC)

    endif

    tonne_to_kgm2s=1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)
    conv = tonne_to_kgm2s !DSRC * EmisDef(eindex_vol)%conv
    do volc_no=1,nvolc 
       i=i_volc(volc_no)
       j=j_volc(volc_no)
       !Find global<->local coordinates for xm2
       if ((i >= gi0).and.(i<=gi1).and.(j>= gj0).and.&
            (j<= gj1))then !on the correct processor
          Volcanoes_found=.true.
          if ( DEBUG_VULC ) write(*,*)'i,j for volcanoe is',i,j
          if ( DEBUG_VULC ) write(*,*)'EMIS_VOLC is',emis_volc(volc_no)
          i_l = i -gi0 +1
          j_l = j- gj0 +1
          if ( DEBUG_VULC ) write(*,*)'Local coord is',i_l,j_l,gi0,gj0
          emis_volc(volc_no) = emis_volc(volc_no)* conv * xm2(i_l,j_l)
       endif
    enddo !volc_no
    !   endif

    !/** broadcast volcano heights
    CALL MPI_BCAST(height_volc,4*NMAX_VOLC,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

  end subroutine VolcGet
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Set_Volc

 !-----------------------------------------------------------------------!
 ! Starts converting emission units from kg/m2/s to.... (hourly)
 !-----------------------------------------------------------------------!

    !**Local variables
    integer            :: k,i,j, i1,i2,j1,j2
    integer            :: iqrc  ! index in emission arrays from EmisGet
    real               :: unit_conv1  

    rcemis_volc0(:) = 0.0
    unit_conv1      = 0.0

    !/** Set volcano
    do volc_no=1,nvolc
       k=height_volc(volc_no)
       i=i_volc(volc_no) +IRUNBEG-1   !NEW
       j=j_volc(volc_no) +JRUNBEG-1   !NEW

       if ( DEBUG_VULC ) &
       write(6,'(a20/4i6/6i6/4i6)')'Volcan: check1 ',  &
       i,j, i_volc(volc_no),j_volc(volc_no),           &
       i_local(i),j_local(j), li0, li1, lj0, lj1,      &
       gi0,gi1,gj0,gj1

       if ( (i_local(i) >= li0) .and. (i_local(i) <= li1 )  .and.  &
            (j_local(j) >= lj0) .and. (j_local(j) <= lj1) ) then 

          unit_conv1 = GRAV* 0.001*AVOG/ &
                      (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))

          rcemis_volc0(volc_no) = emis_volc(volc_no)  &
                                * unit_conv1 / 64.0 !DSRC molwt(QRCVOL)
                                   ! HARD_CODED 64 - fix in future

         if ( DEBUG_VULC ) &
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)

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

   integer i,j,k,i_l,j_l, i1,i2,j1,j2
   real unit_conv2 


  do volc_no=1,nvolc

     k=height_volc(volc_no)
     i=i_volc(volc_no) +IRUNBEG-1   !NEW
     j=j_volc(volc_no) +JRUNBEG-1   !NEW
    ! i=i_volc(volc_no)
    ! j=j_volc(volc_no)

     if ( DEBUG_VULC ) &
     write(6,'(a20/4i6/6i6/4i6)')'Volcan: check2 ', &
     i,j, i_volc(volc_no),j_volc(volc_no),           &
     i_local(i),j_local(j), li0, li1, lj0, lj1,      &
     gi0,gi1,gj0,gj1

       if ( (i_local(i) >= li0) .and. (i_local(i) <= li1 )  .and.  &
            (j_local(j) >= lj0) .and. (j_local(j) <= lj1) ) then 

        i_l = i_local(i) !local i
        j_l = j_local(j) !local j


        if ( DEBUG_VULC ) &
        write(6,'(a30,4i8)')'Volcan: check 3: ',   &
        i_l, j_l, i_volc(volc_no)-gi0+1, j_volc(volc_no)-gj0+1

        unit_conv2 = roa(i_l,j_l,KMAX_BND-k,1) / (ps(i_l,j_l,1)-PT)

        rcemis_volc(volc_no) = rcemis_volc0(volc_no) * unit_conv2

        if ( DEBUG_VULC ) &
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)

     endif
  enddo ! volc_no

  end subroutine Scale_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module Volcanos_ml
