
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!>MODULE  Io_Routines   - small help routines
!! Routines to write out data-arrays in a standardised form. Should help to
!! simplify files for plotting, etc.
!! Very preliminary. Only writedata used so far

   module Io_Routines
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

  use Io_ml, only : IO_LOG 
  implicit none
  private

  public :: writedata
  public :: writetdata
  public :: writeZarray

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine writedata(ofile, headers,coord, data, headerline, odir)

 !/ Writes out data in form which is easy to plot. One
 !/ line of headers (prefix #), then 2-D data array. Needs two calls
!! and reshape command. 

  character(len=*)    :: ofile    !> outout file, will add .txt below
  character(len=*), dimension(:)    :: headers
  real,             dimension(:)    :: coord   !(e.g. z, x)
  real,             dimension(:,:)  :: data
  character(len=*), optional        :: headerline
  character(len=*), optional        :: odir  ! output directory
  integer :: ic, iz, ncols, nrows, io
  character(len=100) :: outfile  ! output directory

  outfile = trim( ofile )// ".txt"
  if( present( odir ) ) outfile = trim(odir)// "/" // trim(ofile)

    nrows = size(data,1)
    ncols = size(data,2)    ! +1 for coord
    print *, "WRITEDATA rows, cols ", nrows, ncols , "FILE ", trim(ofile)

    open(newunit=io,file=trim(ofile) )
    if( present( headerline) ) write(io,"(a)") trim(headerline)
    write(io,"(a,100a12)" )  "#",( trim(headers(ic)), ic=1, ncols+1) ! Header line


    do iz = size(coord), 1, -1  ! Uses z so far, so invert
      write(io,"(f12.2,88es12.3)") coord(iz), ( data(iz, ic), ic=1,ncols)
    end do
    close(io)

end subroutine writedata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> write_tdata  outputs time-varying arrays in form which is easy to plot. One
!! line of headers, then 2-D data array.

subroutine writetdata(ofile, times,coord, data, odir)

  character(len=*)    :: ofile    !> outout file, will add .txt below
  real, intent(in), dimension(:)    :: coord   !(e.g. z, x)
  real, intent(in), dimension(:)    :: times   ! e.g. 0.1, 10.0, 100.
  real, intent(in), dimension(:,:)  :: data
  character(len=*), optional        :: odir  ! output directory

  character(len=100) :: outfile  ! output directory
  real, dimension(size(data,2))  :: odata
  integer :: ic, iz, ncols, nrows, io
  character(len=20), dimension(size(times)+1) :: headers
  character(len=100) :: fmt

  outfile = trim( ofile )// ".txt"
  if( present( odir ) ) outfile = trim(odir)// "/" // trim(ofile)
  print *, "OUTFILE BBB", trim(outfile)
!stop 'BBB'

    nrows = size(data,1)
    ncols = size(data,2)    ! +1 for coord

    !open(newunit=io,file="Results_" // trim(label) // ".txt")
    open(newunit=io,file=trim(ofile) ) ! // ".txt" )

   ! Headers have times
    headers(1) = "0"   !! "t:"
    do ic = 1, size(times)
       fmt="(f10.3)"
       if( times( ic) > 1.0 ) fmt="(f10.2)"
       if( times( ic) > 10.0 ) fmt="(f10.1)"
       write(headers(ic+1),fmt=fmt) times( ic)  ! convert to txt
    end do
    write(io,"(100a10)" )  ( trim(headers(ic)), ic=1, ncols+1) ! Header line

    !print *, "WRITEDATA rows, cols ", nrows, ncols 

    !INVERSE do iz = size(coord), 1, -1  ! Uses z so far, so invert
    do iz = 1, size(coord)
       odata = min( data(iz,:) , 9.99e29 ) !! Looking for v. strange values!
     !! gfortran printed e.g. 1.00e-2345 as *****. Reset, and mark with -ve:
      !! where(odata < 1.0e-99 .and. odata > 0) odata = -0.0
      where(odata < tiny(1.0) .and. odata > 0) odata = -0.0
      !write(io,"(f10.2,88es10.2)") coord(iz), ( odata(ic), ic=1,ncols)
      write(io,"(f10.2,88es12.4)") coord(iz), ( odata(ic), ic=1,ncols)
      !write(io,*) coord(iz), ( data(iz, ic), ic=1,ncols)
    end do
    close(io)

end subroutine writetdata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine writeZarray(io,label, coord, tcoord, data)

 !/ Writes out data in form which is easy to plot. One
 !/ line of z-coords, then 2-D data array - 1st column should be time.

  integer, intent(in) :: io
  character(len=*)    :: label
  real,             dimension(:)    :: coord   !(e.g. z, x)
  real,             dimension(:)    :: tcoord  ! time-coord
  real,             dimension(:,:)  :: data
  integer :: ir, iz, ncols, nrows

    nrows = size(data,2)
    ncols = size(coord)    ! +1 for coord

    open(io,file="ZZarray_" // trim(label) // ".txt")
  ! write z coords in row header. Start with 0 - will match time.coord
    write(io,"(a10,99f7.1)" )  " 0 " , ( coord(iz), iz=1, ncols ) ! Header line

    print *, "WRITEDATA rows, cols ", nrows, ncols 

    do ir = nrows, 1, -1
      write(io,"(f10.2,88f7.1)") tcoord(ir), ( data(iz, ir), iz=1,ncols)
      !print   *, "PRINT: ", coord(iz), ( data(iz, ic), ic=1,ncols)
    end do
    
    close(io)

end subroutine writeZarray

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
                        end module Io_Routines
! __________________________________________________________________________ !
