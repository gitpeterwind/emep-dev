                         module Par_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define all parameters and variables used in domain definition and
!  the parallel data decomposition.  
!
! uni - now includes
! eulpar.inc
! eulnx.inc
!----------------------------------------------------------------------------
!  $Id: Par_ml.pat,v 1.3 2002-09-13 08:49:00 mifads Exp $
!  Erik Berge, DNMI    Roar Skaalin, SINTEF Industrial Mathematics
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!     Parameter statements for definition of integration area
!     and maximum number of gridpoints per processor
!RESTRI
!  we try to run on a smaller domain with the same input
!  for this reason we now define additional parameters:
!  numbers of points in the 'larger' array iilardom,jjlardom
!  and coordinates of the origin of the smaller domain with
!  respect to the larger domain ismbeg,jsmbeg
!  one can run the large domain by setting  
!  gimax, gjmax to iilardom,jjlardom
!  and  ismbeg,jsmbeg   to 1
!  also we have now to distinguish mfsize for input and output: 
!  mfsizeinp,mfsizeout!!!
!RESTRI
!
!RESTRI + uni

implicit none
private

  integer, public, parameter ::  &
    IILARDOM    =    170         &
  , JJLARDOM    =    133         &
  , ISMBEG      =   domainx0     &
  , JSMBEG      =   domainy0     &
  , GIMAX       =   domaindx     & ! Number of global points in longitude
  , GJMAX       =   domaindy     & ! Number of global points in latitude
  , NPROCX      =   nprocx       & ! Actual number of processors in longitude
  , NPROCY      =   nprocy       & ! Actual number of processors in latitude
!RESTRI + uni
  , MAXLIMAX   = (GIMAX+NPROCX-1)/NPROCX   & ! Maximum number of local points in longitude
  , MAXLJMAX   = (GJMAX+NPROCY-1)/NPROCY   & ! Maximum number of local points in latitude
!RESTRI
  , MFSIZEINP = IILARDOM*JJLARDOM   & ! Maximum field size for input
  , MFSIZEOUT = GIMAX*GJMAX           ! Maximum field size for output
!RESTRI
!
 integer , public, parameter :: &
      NPROC = NPROCX*NPROCY     ! Actual total number of processors
!

!
!     Parameter statements for the parameters used to access the tabell 
!     of neighbor processors (neighbor)
!
 integer , public, parameter :: &
     NORTH = 1   &   ! Neighbor to the north
   , SOUTH = 2   &   ! Neighbor to the south
   , EAST  = 3   &   ! Neighbor to the east
   , WEST  = 4   &   ! Neighbor to the west
   , NOPROC = -1     ! Value in neighbor when there is no neighbor in the
                          ! actual direction
!
!
!     Variables for actual number of local points, to be computed
!
 integer , public, save :: &
     limax       & ! Actual number of local points in longitude
   , ljmax         ! Actual number of local points in latitude

!
!     Variables for global address of the start and end points 
!     on each processor
!
 integer, public, save :: & 
     gi0        &  ! Global address of longitude start point
   , gi1        &  ! Global address of longitude end point
   , gj0        &  ! Global address of latitute start point
   , gj1           ! Global address of latitude end point

!
!     Variables used as loop indexes on each processor. The values are
!     derived from limax and ljmax and the position of the domain belonging
!     to this processor. See PARINIT for assignment of values.
!
 integer, public, save :: &
     li0  & ! First local index in longitude when outer boundary is excluded
   , li1  & ! Last local index in longitude when outer boundary is excluded
   , lj0  & ! First local index in latitude when outer boundary is excluded
   , lj1    ! Last local index in latitude when outer boundary is excluded
!
!
!     Variables for address of this processor
!
 integer , public, save :: &
   me   &  ! Address of processer, numbering starts at 0 in south-west corner of ground level
 , mex  &  ! Longitude address of processor, numbering  starts at 0 on the westmost boundary
 , mey     ! Latitude address of processor, numbering starts at 0 on the southmost boundary

!
!     Variable for the table of local neighbor
!

 integer, public, save, dimension(4) ::  neighbor

!
!     Tables of actual number of points and start and end points for
!     all processors
!
 integer, public, save, dimension(0:NPROC-1) ::  &
      tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
!
!
!     All variables *which were* stored in common blocks PARTOPO (single variables 
!     and local tables) and TPARTOPO (tables common to all processors)
!
! integer, public, save ::                   &
!            limax,  li0,  li1,  gi0,  gi1   &
!         ,  ljmax,  lj0,  lj1,  gj0,  gj1   &
!        ,  me,     mex,    mey
! defined above:         ,  neighbor

! integer, public, save ::       & 
!      tlimax,  tgi0,  tgi1      &
!    , tljmax,  tgj0,  tgj1


!--- eulnx.inc follows:

!------------------------------------------------------------------------------
! Define parameters used in the communication
!------------------------------------------------------------------------------
!
!     Code for broadcasting information to all nodes
!
    integer, public, parameter ::  NX_BROADCAST = -1
!
!     The different messages used in the bott version of airpol
!
    integer, public, parameter :: &
               MSG_INIT0 = 10     &
              ,MSG_INIT1 = 11     &
              ,MSG_INIT2 = 12     &
              ,MSG_INIT3 = 13     &
              ,MSG_INIT4 = 14     &
              ,MSG_INIT5 = 15     &
              ,MSG_INIT6 = 16     &
              ,MSG_INIT7 = 17     &
              ,MSG_INIT8 = 18     &
              ,MSG_INIT9 = 19     &
              ,MSG_NORTH1 = 21     &
              ,MSG_NORTH2 = 22     &
              ,MSG_EAST1 = 31     &
              ,MSG_EAST2 = 32     &
              ,MSG_SOUTH1 = 41     &
              ,MSG_SOUTH2 = 42     &
              ,MSG_WEST1 = 51     &
              ,MSG_WEST2 = 52     &
              ,MSG_TOPO1 = 61     &
              ,MSG_TOPO2 = 62  
    integer, public, parameter :: & 
               MSG_MAIN1 = 71     &
              ,MSG_MAIN2 = 72     &
              ,MSG_READ1 = 81     &
              ,MSG_READ2 = 82     &
              ,MSG_READ3 = 83     &
              ,MSG_READ4 = 84     &
              ,MSG_READ5 = 85     &
              ,MSG_READ7 = 87     &
              ,MSG_FIELD1 = 91     &
              ,MSG_MET1 = 101     &
              ,MSG_MET2 = 102     &
              ,MSG_MET3 = 103     &
              ,MSG_MET4 = 104     &
              ,MSG_MET5 = 105     &
              ,MSG_MET6 = 106     &
              ,MSG_PARI = 107
!-- end of eulnx.inc

 public :: parinit

 contains

	subroutine parinit(min_grids)   !u4 - argument added to avoid "use"

!	use PAR_ML   , only: GIMAX,GJMAX,NPROCX,NPROCY,NPROC   &
!     			,li0,li1,limax,lj0,lj1,ljmax           &
!     			,gi0,gi1,gj0,gj1		       &
!     			,me,mex,mey			       &
!     			,tgi0,tgi1,tlimax,tgj0,tgj1,tljmax     &
!     			,neighbor,EAST,WEST,SOUTH,NORTH,NOPROC
!u4 	use Advection_ml , only: MIN_ADVGRIDS

	implicit none
        integer, intent(in) :: min_grids  ! u4
	integer i, j, ime, imex, imey, rest
           
!
!
!     Find the x-, y-, and z-addresses of the domain assigned to the
!     processor
!
!     pw 09.01.2002: Check if the subdomain is large enough
!     pw 26.03.2002: Bug in formula for tljmax, tgj0, tgj1  
!
	mey = me/NPROCX
	mex = me - mey*NPROCX
!
!
!     Find the number of grid points in each direction for this processor.
!     We first try to divide the total number equally among the 
!     processors. Then the rest is distributed one by one to first processor 
!     in each direction. Here we also set the global address of the start
!     and end point in each direction.
!
!     x-direction (longitude)
!
	limax = GIMAX/NPROCX
	rest = GIMAX - limax*NPROCX
	gi0 = mex*limax + 1
	if(rest>0)then
	  if(mex.eq.NPROCX-1)then
	    limax = limax+1
	    gi0 = gi0+rest-1
	  elseif (mex < rest-1) then
	    limax = limax + 1
	    gi0 = gi0 + mex
	  else
	    gi0 = gi0 + rest-1
	  endif
	endif
	gi1 = gi0 + limax - 1
!
!     y-direction (latitude)
!
	ljmax = GJMAX/NPROCY
	rest = GJMAX - ljmax*NPROCY
	gj0 = mey*ljmax + 1
	if(rest>0)then
	  if(mey.eq.NPROCY-1)then
	    ljmax = ljmax + 1
	    gj0 = gj0 + rest-1
	  elseif (mey < rest-1) then
	    ljmax = ljmax + 1
	    gj0 = gj0 + mey
	  else
	    gj0 = gj0 + rest-1
	  endif
	endif
	gj1 = gj0 + ljmax - 1
!
!     Find the neighbors of this processor.
!     NB! Assumes non-cyclic boundary conditions
!
	if (mey > 0) then
	  neighbor(SOUTH) = me-NPROCX
	  lj0 = 1
	else
	  neighbor(SOUTH) = NOPROC
	  lj0 = 2
	endif
	if (mey < NPROCY-1) then
	  neighbor(NORTH) = me+NPROCX
	  lj1 = ljmax
	else
	  neighbor(NORTH) = NOPROC
	  lj1 = ljmax - 1
	endif
	if (mex > 0) then
	  neighbor(WEST) = me-1
	  li0 = 1
	else
	  neighbor(WEST) = NOPROC
	  li0 = 2
	endif
	if (mex < NPROCX-1) then
	  neighbor(EAST) = me+1
	  li1 = limax
	else
	  neighbor(EAST) = NOPROC
	  li1 = limax - 1
	endif
!
!
!     Initialize the tables containing number of gridpoints and
!     addresses of start and endpoint in all directions, for all
!     processors. This is a repetition of the computations above,
!     but now for all processors.
!
	do ime = 0, NPROC-1
!
	  imey = ime/NPROCX
	  imex = ime - imey*NPROCX
!
!     x-direction (longitude)
!
	  tlimax(ime) = GIMAX/NPROCX
	  rest = GIMAX - tlimax(ime)*NPROCX
	  tgi0(ime) = imex*tlimax(ime) + 1
	  if(rest>0)then
	    if (imex .eq. NPROCX-1) then
	      tlimax(ime) = tlimax(ime) + 1
	      tgi0(ime) = tgi0(ime) + rest-1
	    elseif (imex < rest-1) then
	      tlimax(ime) = tlimax(ime) + 1
	      tgi0(ime) = tgi0(ime) + imex
	    else
	      tgi0(ime) = tgi0(ime) + rest-1
	    endif
	  endif
	  tgi1(ime) = tgi0(ime) + tlimax(ime) - 1
!
!     y-direction (latitude)
!
	  tljmax(ime) = GJMAX/NPROCY
	  rest = GJMAX - tljmax(ime)*NPROCY
	  tgj0(ime) = imey*tljmax(ime) + 1
	  if(rest > 0)then
	    if (imey .eq. NPROCY-1) then
	      tljmax(ime) = tljmax(ime) + 1
	      tgj0(ime) = tgj0(ime) + rest-1
! pw u3	    elseif (imey .le. rest-1) then
	    elseif (imey < rest-1) then
	      tljmax(ime) = tljmax(ime) + 1
	      tgj0(ime) = tgj0(ime) + imey
	    else
	      tgj0(ime) = tgj0(ime) + rest-1
	    endif
	  endif
	  tgj1(ime) = tgj0(ime) + tljmax(ime) - 1
!
	enddo
!
!	write(6,111) me, mex, mey, neighbor, gi0, gj0, limax, ljmax        
! 111	format(11i7)

!       pw The size of the grid cannot be too small.       

	!u4 if(limax < MIN_ADVGRIDS .or. ljmax < MIN_ADVGRIDS )then
	if(limax < min_grids .or. ljmax < min_grids )then

	   write(*,*)'subdomain too small!'
	   write(*,*)'lXmax must be at least min_grids = ', &
     	                                        min_grids
	   write(*,*)'me,li0,li1,limax ',me,li0,li1,limax
	   write(*,*)'me,lj0,lj1,ljmax ',me,lj0,lj1,ljmax
	   call gc_abort(me,NPROC,'SUBDOMAIN TO SMALL!!!')
	endif

	return
	end subroutine parinit

end module Par_ml
