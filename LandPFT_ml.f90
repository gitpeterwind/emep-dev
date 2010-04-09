! <LandPFT_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module LandPFT_ml

use CheckStop_ml,   only: CheckStop
use GridValues_ml,  only: debug_proc, debug_li, debug_lj
use ModelConstants_ml,  only : DEBUG_LANDPFTS, MasterProc
use NetCDF_ml, only: ReadField_CDF !dsLPJ
use Par_ml,         only: MAXLIMAX, MAXLJMAX
use SmallUtils_ml,  only: find_index, NOT_FOUND, WriteArray

implicit none
private


!/- subroutines:

  public :: MapPFT_LAI
  public :: MapPFT_BVOC


 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 character(len=80), private :: errmsg

 logical, public, parameter :: PFTS_USED = .false.
 real, public, allocatable :: pft_lai(:,:,:) 
 real, public, allocatable :: pft_bvoc(:,:,:,:) 

 ! PFTs available from smoothed LPJ fields

  integer, public, parameter :: N_PFTS = 6
  character(len=5),public, parameter, dimension(N_PFTS) :: PFT_CODES = &
        (/ "CF   ", "DF   ", "NF   ", "BF   ", "C3PFT", "C4PFT" /) 

   ! Variables available:

    character(len=5),public, parameter :: LAI_VAR =  "LAIv_"
    character(len=5),public, parameter, dimension(2) :: BVOC_VAR = &
        (/ "Eiso_" , "Emt_ " /)

   !skip (/ "Normed_LAIv", "LAIv_      ", "Emt_       ", "Eiso_      " /)

contains

 !==========================================================================
 subroutine MapPFT_LAI(month)

!.....................................................................
!**    DESCRIPTION:

!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!    The LPJ data have been merged into 4 EMEP forest classes and two
!    other veg, for either C3 or C4 vegetation.
!    Normed_LAIv is relative LAI, with max value 1.0


    integer, intent(in) :: month

    real    :: lpj(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    logical :: my_first_call = .true.
    integer ::  n, pft, ivar, iLC, alloc_err, iiLC
    character(len=20) :: varname

     if ( my_first_call ) then
         allocate ( pft_lai(MAXLIMAX,MAXLJMAX,N_PFTS) )
         my_first_call = .false.
     end if
         
    ! Get LAI data:

     do pft =1, N_PFTS
           varname = trim(LAI_VAR) // trim(PFT_CODES(pft)) 

           call ReadField_CDF('GLOBAL_LAInBVOC.nc',varname,&
              lpj,month,interpol='zero_order',needed=.true.,debug_flag=.true.)

           pft_lai(:,:,pft ) = lpj(:,:)

     end do ! pft

  end subroutine MapPFT_LAI

 !==========================================================================

 subroutine MapPFT_BVOC(month,bvoc_wanted)

!.....................................................................
!**    DESCRIPTION:

!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!    The LPJ data have been merged into 4 EMEP forest classes and two
!    other veg, for either C3 or C4 vegetation.
!    Normed_LAIv is relative LAI, with max value 1.0


    integer, intent(in) :: month
    character(len=*), dimension(:), intent(in) :: bvoc_wanted   ! e.g. "C5H8"

    real    :: lpj(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    logical :: my_first_call = .true.
    integer ::  n, pft, ivar
    character(len=20) :: varname

    ! Ebvoc already includes monthly LAI changes - might be wrong?
    
     if ( my_first_call ) then
         allocate ( pft_bvoc(MAXLIMAX,MAXLJMAX,N_PFTS,size(bvoc_wanted)) )
         my_first_call = .false.
     end if
         
    ! Get BVOC data. Code assumes that we want isoprene first, then
    !    apinene if provided. Multiple terpenes not considered yet,
    !    but we have just Emt anyway.


     do pft =1, N_PFTS
       do ivar =1, size( bvoc_wanted ) 
           varname = trim(BVOC_VAR(ivar)) // trim(PFT_CODES(pft)) 

           call ReadField_CDF('GLOBAL_LAInBVOC.nc',varname,&
              lpj,month,interpol='zero_order',needed=.true.,debug_flag=.true.)

           pft_bvoc(:,:,pft, ivar ) = lpj(:,:)

       end do ! ivar
     end do ! pft


  end subroutine MapPFT_BVOC

 !=======================================================================
end module LandPFT_ml
