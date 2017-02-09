!
! NAME
!   file_nc  -   tools to write NetCDF files
!
! USAGE
!
!   use file_nc
!
!   !
!   ! dump a single data set into a NetCDF file;
!   ! input arrray 'dat' is integer or real with rank up to 7;
!   ! provide data set name and dimension names;
!   ! return status non-zero in case of any error;
!   ! optional character key 'typ' to convert to other type:
!   !   'int', 'integer', 'integer(4)'
!   !   'real', 'float', 'real(4)'
!   !   'double', 'real(8)'
!   !
!   call nc_dump( 'data.nc', dat, 'dat', (/'x','y','z'/), status [,typ='real(4)'] )
!
! AUTHOR
!   Arjo Segers, TNO, The Netherlands
!

module file_nc

  use file_nc_i4
  use file_nc_r4
  use file_nc_r8

  implicit none

 
  ! --- in/out ------------------------------------
  
  public

end module file_nc
