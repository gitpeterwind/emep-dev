module Sites_ml
! -----------------------------------------------------------------------
! Contains subroutines to read in list of measurement stations and/or 
! radiosonde locations. 
!
! Developed from sitesout.f. Converted to  F90, and written so that sites
! are now specified in the input file "sites.dat" and "sondes.dat"
!
! Feb-Mar/2001 Cleaning up, removal of stop_test,MAP_ADV,  hf/ds
! 4-5/2001 Substantial changes to harmonise treatment of sites and sondes.
!     and to enable flexible choice of outputs in connection with
!     the My_Outputs_ml.
! 26/4/01 - Sonde read-in added to sitesdef.
! 24/4/01 - explicit description data written to sites.out
! -----------------------------------------------------------------------

use My_Outputs_ml, only : &   ! for sitesout
        NSITES_MAX, &
          NADV_SITE, NSHL_SITE, NXTRA_SITE, &
          SITE_ADV, SITE_SHL, SITE_XTRA, &
          FREQ_SITE, &
        NSONDES_MAX, &
          NADV_SONDE, NSHL_SONDE, NXTRA_SONDE, &
          SONDE_ADV, SONDE_SHL, SONDE_XTRA, &
          FREQ_SONDE
use My_Derived_ml, only : d_2d, IOU_INST, &
          D2_VG_REF, D2_VG_1M, D2_VG_STO, D2_FX_REF, D2_FX_STO

use Par_ml , only : ISMBEG,JSMBEG,GIMAX,GJMAX,  &
              GI0,GI1,GJ0,GJ1,me,NPROC,MAXLIMAX,MAXLJMAX
use Io_ml  , only : check_file,open_file,ios,fexist,IO_SITES,IO_SONDES
use GenSpec_adv_ml   , only : NSPEC_ADV
use GenSpec_shl_ml   , only : NSPEC_SHL
use GenChemicals_ml  , only : species                    ! for species names
use ModelConstants_ml, only : NMET,PPBINV,PPTINV,KMAX_MID, current_date,KMAX_BND
implicit none
private                             ! stops variables being accessed outside


!/** subroutines made available ** /

public :: sitesdef               ! Calls Init_sites for sites and sondes
public :: siteswrt_surf          ! Gets site  data ready for siteswrt_out
public :: siteswrt_sondes        ! Gets sonde data ready for siteswrt_out
private :: Init_sites            ! reads locations, species
private :: set_species           ! Sets species/variable names for output
private :: siteswrt_out          ! Collects output from all nodes and prints


!/** some variables used in following subroutines **/

integer, private, save :: nglobal_sites, nlocal_sites
integer, private, save :: nglobal_sondes, nlocal_sondes

! /**  site_gindex stores the global index n asociated with each 
!      processor and local site   **/

integer, private, save, dimension (0:NPROC-1,NSITES_MAX) :: site_gindex
integer, private, save, dimension (0:NPROC-1,NSONDES_MAX) :: sonde_gindex

integer, private, save, dimension (NSITES_MAX) ::  &
         site_gx, site_gy   & ! global coordinates
       , site_x, site_y     & ! local coordinates
       , site_n               ! number in global
integer, private, save, dimension (NSONDES_MAX) ::  &
         sonde_gx, sonde_gy   & ! global coordinates
       , sonde_x, sonde_y     & ! local coordinates
       , sonde_n               ! number in global

! Values from My_Outputs_ml gives ... => 

integer, public, parameter :: &
     NOUT_SITE =            NADV_SITE  + NSHL_SITE + NXTRA_SITE  &! Total No.
    ,NOUT_SONDE = KMAX_MID * ( NADV_SONDE + NSHL_SONDE+ NXTRA_SONDE )

character(len=50), private, save, dimension (NSITES_MAX) ::  site_name
character(len=50), private, save, dimension (NSONDES_MAX)::  sonde_name
character(len=20), private, save, &
              dimension (NADV_SITE+NSHL_SITE+NXTRA_SITE) :: site_species
character(len=20), private, save, &
           dimension (NADV_SONDE+NSHL_SONDE+NXTRA_SONDE) :: sonde_species

character(len=40), private  :: errmsg   ! Message text
integer, private,  save  :: &   !u3 - global versions added
       gibegpos,giendpos,gjbegpos,gjendpos &!/** domain for global coords **/ 
      ,ibegpos,iendpos,jbegpos,jendpos      !/** domain for local node    **/
integer, private  :: d, info            ! processor index and gc_send info
integer, private  :: i, n, nloc, ioerr  ! general integers

 !-- Debugging parameter:
 logical, private, parameter :: MY_DEBUG = .false.

contains

!==================================================================== >
   subroutine sitesdef()

   ! -------------------------------------------------------------------------
   ! reads in sites.dat and sondes.dat (if present), assigns sites to
   ! local domains, and collects lists of sites/species/variables for output.
   ! -------------------------------------------------------------------------

     call Init_sites("sites",IO_SITES,NSITES_MAX, &
           nglobal_sites,nlocal_sites, &
           site_gindex, site_gx, site_gy, site_x, site_y, site_n, &
           site_name)

     call Init_sites("sondes",IO_SONDES,NSONDES_MAX, &
           nglobal_sondes,nlocal_sondes, &
           sonde_gindex, sonde_gx, sonde_gy, sonde_x, sonde_y, sonde_n, &
           sonde_name)

     call set_species(SITE_ADV,SITE_SHL,SITE_XTRA,site_species)
     call set_species(SONDE_ADV,SONDE_SHL,SONDE_XTRA,sonde_species)

     if ( MY_DEBUG ) then
        write(6,*) "sitesdef After nlocal ", nlocal_sites, " on me ", me 
        do i = 1, nlocal_sites
          write(6,*) "sitesdef After set_species x,y ", &
                        site_x(i), site_y(i), " on me ", me
        end do
     end if ! DEBUG

   end subroutine sitesdef
!==================================================================== >
   subroutine set_species(adv,shl,xtra,s)

   ! -------------------------------------------------------------------------
   ! Makes a character array "s" containg the names of the species or
   ! meteorological parameters to be output. Called for sites and sondes.
   ! -------------------------------------------------------------------------
   
      integer, intent(in), dimension(:) :: adv, shl ! Arrays of indices wanted
      character(len=*), intent(in), &
                             dimension(:) :: xtra  !Names of extra params
      character(len=*),intent(out), dimension(:) :: s

      integer :: nadv, nshl, n2, nout     ! local sizes
      nadv = size(adv) 
      nshl = size(shl) 
      n2 = nadv + nshl
      nout = size(s)             ! Size of array to be returned

       s(1:nadv)    = species( NSPEC_SHL + adv(:) )%name
       s(nadv+1:n2) = species( shl(:) )%name
       s(n2+1:nout) = xtra(:)

  end  subroutine set_species

!==================================================================== >
   subroutine Init_sites(fname,io_num,NMAX, &
           nglobal,nlocal, &
           s_gindex, s_gx, s_gy, s_x, s_y, s_n, s_name)
   ! -------------------------------------------------------------------------
   ! Reads the file "sites.dat" and "sondes.dat" to get coordinates of 
   ! surface measurement stations or locations where vertical profiles
   ! or extra output are required. (These files may be empty, but this is
   ! not recommended - the sites data provide good diagnostics).
   !RESTRI
   !  define, whether a certain output site belongs to the given processor
   !  and assign the local coordinates
   !RESTRI
   !
   ! -------------------------------------------------------------------------
   ! NB. global below refers to all nodes
   !     local  below refers to the local node

   character(len=*), intent(in) :: fname
   integer,          intent(in) :: io_num
   integer,          intent(in) :: NMAX       ! Max no. sites
   integer, intent(out)  :: nglobal, nlocal   ! No. sites 
   integer, intent(out), dimension (0:,:) :: s_gindex  ! index, starts at me=0
   integer, intent(out), dimension (:) ::  &   
                              s_gx, s_gy   & ! global coordinates
                            , s_x, s_y     & ! local coordinates
                            , s_n            ! number in global
   character(len=*), intent(out), dimension (:) ::  s_name
 
   !-- Local:
   integer,  dimension (NMAX) :: s_n_recv  ! number in global

   integer           :: nin      ! loop index
   real              :: x, y     ! coordinates read in
   character(len=20) :: s        ! Name of site read in
   character(len=30) :: comment  ! comment on site location
   character(len=40) :: infile

    infile  = fname // ".dat"

    call check_file(infile,fexist,needed=.false.)

     if ( .not. fexist ) return

    !/** intialise RESTRI domain
     n       = 0                   ! No. sites found within domain
     gibegpos = ISMBEG
     giendpos = GIMAX+ISMBEG-1
     gjbegpos = JSMBEG
     gjendpos = GJMAX+JSMBEG-1

   ios = 0   ! zero indicates no errors
   errmsg = "ios error" // infile

   if(me == 0) call open_file(io_num,"r",infile,needed=.true.)

   if ( NMAX /= size(s_name) ) then   !-- consistency check
          call gc_abort(me,NPROC,"sitesdef NMAX problem")
   end if

   

   SITEREAD: if(me == 0) then   
     SITELOOP:  do nin = 1, NMAX

         read (unit=io_num,fmt=*,iostat=ioerr) s, x, y

           if ( ioerr < 0 ) then
             write(6,*) "sitesdef : end of file after ", nin-1,  infile
             exit SITELOOP
           end if ! ioerr

         !RESTRI:  check if the output site belongs to the restricted domain 
         !         at all, if not print a warning message

         if(   gibegpos > x .or. giendpos < x .or.  & 
               gjbegpos > y .or. gjendpos < y )  then
               write(6,*) "sitesdef: ", s, x,y, " outside computational domain"
         else
              n = n + 1
              s_gx(n)   = x  
              s_gy(n)   = y  
              if ( x == gibegpos .or. x == giendpos .or. &
                   y == gjbegpos .or. y == gjendpos ) then

                    comment = " WARNING - domain boundary!!"
              else
                    comment = " ok - inside domain         "
              end if
              s_name(n)  = s // comment           

         endif


  
     end do SITELOOP 

     nglobal       = n 

    !/ NSITES/SONDES_MAX must be _greater_ than the number used, for safety
    !  (we could in fact dimesnion as NMAX-1 but I can't be bothered..)

     if ( n >= NMAX ) call gc_abort(me,NPROC,"increase NGLOBAL_SITES_MAX")

     close(unit=io_num)
   end if SITEREAD


   ! - send global coordinates to processors
   call gc_ibcast(680,1,0,NPROC,info,nglobal)
   call gc_ibcast(681,nglobal,0,NPROC,info,s_gx)
   call gc_ibcast(682,nglobal,0,NPROC,info,s_gy)

!/**   define local coordinates of first and last element of the
!      arrays with respect to the larger domain  **/

    ibegpos = gi0+ISMBEG-1
    iendpos = gi1+ISMBEG-1
    jbegpos = gj0+JSMBEG-1
    jendpos = gj1+JSMBEG-1


    nlocal   = 0

    do n = 1, nglobal

      if(ibegpos <= s_gx(n) .and.  iendpos >= s_gx(n) .and.   &
         jbegpos <= s_gy(n) .and.  jendpos >= s_gy(n) ) then

        nlocal      = nlocal + 1
        s_x(nlocal) = s_gx(n)-ibegpos+1
        s_y(nlocal) = s_gy(n)-jbegpos+1
        s_n(nlocal) = n

        if ( MY_DEBUG ) then
           write(6,*) "sitesdef Site on me : ", me, " No. ", nlocal, &
                s_gx(n), s_gy(n) , " =>  ", &
                s_x(nlocal), s_y(nlocal)
        end if

      endif
    end do ! nglobal

    !-- inform me=0 of local array indices:
    if (MY_DEBUG) write(6,*) "sitesdef ", fname, " before gc NLOCAL_SITES", me, nlocal

    if ( me .ne. 0 ) then
       if(MY_DEBUG)  write(6,*) "sitesdef ", fname, " send gc NLOCAL_SITES", me, nlocal
       call gc_isend(333,1,0,info,nloc,nlocal)
       if(nlocal > 0)   &
         call gc_isend(334,nlocal,0,info,s_n_recv,s_n)

    else
       if(MY_DEBUG)  write(6,*) "sitesdef for me =0 OCAL_SITES", me, nlocal
       do n = 1, nlocal
         s_gindex(me,n) = s_n(n)
       end do
         
       do d = 1, NPROC-1
          call gc_irecv(333,1,d,info,nloc,nloc)
          if(nloc > 0)   &
         call gc_irecv(334,nloc,d,info,s_n_recv,s_n)
          if (MY_DEBUG) write(6,*) "sitesdef: recv d ", fname, d,  &
                   " zzzz nloc : ", nloc, " zzzz me0 nlocal", nlocal
          do n = 1, nloc
             s_gindex(d,n) = s_n_recv(n)
             if ( MY_DEBUG ) write(6,*) "sitesdef: for d =", fname, d, &
             " nloc = ", nloc, " n: ",  n,  " gives nglob ", s_gindex(d,n)
          end do ! n
       end do ! d
    end if ! me

    if ( MY_DEBUG ) write(6,*) 'sitesdef DS on me', me, ' = ', nlocal


end subroutine Init_sites

!==================================================================== >

  subroutine siteswrt_surf(xn_adv,cfac,xn_shl)
  ! ---------------------------------------------------------------------
  ! writes out just simple concentrations for now....
  ! will be improved later to allow choice of output parameter
  ! should look at chemint also - seems similar for somethings
  ! ---------------------------------------------------------------------
   !u7.4vg use Met_ml, only : temp2m, th   ! Output with concentrations
   use Met_ml, only : t2, th, pzpbl   ! Output with concentrations

  ! -- arguments
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_adv
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX), intent(in)       :: cfac
  real, dimension(NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_shl


  ! Local
  integer :: nglob, nloc, ix, iy, ispec   !/** Site indices           **/
  integer :: nn                           !/** species index          **/
  logical, save :: my_first_call = .true. !/** for debugging          **/

   real,dimension(NOUT_SITE,NSITES_MAX) :: out    !/** for output, local node**/

     if ( MY_DEBUG ) then 
        print *, "sitesdef Into surf  nlocal ", nlocal_sites, " on me ", me 
        do i = 1, nlocal_sites
           print *, "sitesdef Into surf  x,y ",site_x(i),site_y(i)," me ", me
        end do
 
        if ( me == 0 ) then
           print *,  "======= site_gindex ======== sitesdef ============"
           do n = 1, nglobal_sites 
              print "(a12,i4, 2x, 8i4)", "sitesdef ", n, &
                       (site_gindex(d,n),d=0, NPROC-1)  
           end do
           print *,  "======= site_end    ======== sitesdef ============"
        end if ! me = 0
     end if  ! MY_DEBUG

  !/** assign local data to out **/

  do i = 1, nlocal_sites
     ix = site_x(i)
     iy = site_y(i)

     do ispec = 1, NADV_SITE
        out(ispec,i)  = xn_adv( SITE_ADV(ispec) ,ix,iy,KMAX_MID ) * &
                            cfac( SITE_ADV(ispec),ix,iy) * PPBINV

     end do
     my_first_call = .false.

     do ispec = 1, NSHL_SITE
        out(NADV_SITE+ispec,i)  = xn_shl( SITE_SHL(ispec) ,ix,iy,KMAX_MID )
     end do

    !/** then print out XTRA stuff, usually the temmp
    !    or pressure
    !SITE_XTRA=  (/ "th  ", "hmix", "Vg_ref", "Vg_1m", "Vg_sto", "Flux_ref", "Flux_sto" /)

    do ispec = 1, NXTRA_SITE
       nn = NADV_SITE + NSHL_SITE + ispec
       select case ( SITE_XTRA(ispec) )
       case ( "T2" ) 
          out(nn,i)   = t2(ix,iy) - 273.15 
       case ( "th" ) 
          out(nn,i)   = th(ix,iy,KMAX_MID,1)
       case ( "hmix" ) 
          out(nn,i)   = pzpbl(ix,iy)
       case ( "Vg_ref" ) 
          out(nn,i)   = d_2d(D2_VG_REF,ix,iy,IOU_INST)
       case ( "Vg_1m" ) 
          out(nn,i)   = d_2d(D2_VG_1M ,ix,iy,IOU_INST)
       case ( "Vg_sto" ) 
          out(nn,i)   = d_2d(D2_VG_STO ,ix,iy,IOU_INST)
       case ( "Flux_ref" ) 
          out(nn,i)   = d_2d(D2_FX_REF ,ix,iy,IOU_INST)
       case ( "Flux_sto" ) 
          out(nn,i)   = d_2d(D2_FX_STO ,ix,iy,IOU_INST)
       end select 
    end do
  end do

  !/** collect data into gout on me=0 t **/

   call siteswrt_out("sites",IO_SITES,NOUT_SITE,NSITES_MAX, FREQ_SITE, &
                      nglobal_sites,nlocal_sites, &
                     site_gindex,site_name,site_gx,site_gy,site_species,out)

end subroutine siteswrt_surf


!==================================================================== >

  subroutine siteswrt_sondes(xn_adv,xn_shl)
   use Met_ml, only : z_bnd,z_mid,roa,th,xksig   ! Output with concentrations
  ! -------------------------------------------------------------------
  !  Writes vertical concentration  data to files.
  !  IO_SONDES is set in io_ml to be 30
  ! -------------------------------------------------------------------

!  real, dimension(:,:,:),      intent(in) :: z_bnd, z_mid
!  real, dimension(:,:,:,:),    intent(in) :: roa

  real, dimension(:,:,:,:),    intent(in) ::  xn_adv
  real, dimension(:,:,:,:),    intent(in) ::  xn_shl
 
  ! --  Output variables - none

  ! --  Local variables
  integer :: n, i, k,  ix, iy, nn, ispec   ! Site and chem indices


  real, dimension(NOUT_SONDE,NSONDES_MAX) :: out
 
  ! --- ////////////// code ////////////////////////////////// ---

      !MY_DEBUG print *, "Into sondesdef", me, nlocal_sondes

     !** Consistency check 

      do ispec = 1, NXTRA_SONDE
            select case ( SONDE_XTRA(ispec) )
            case ( "z_mid", "xksig ", "th   " )
               errmsg = "ok"
            case default
              errmsg = "ERROR: SONDE_XTRA:" // SONDE_XTRA(ispec) 
!hf u2              call stop_test(.false.,me,NPROC,99,errmsg)
              call gc_abort(me,NPROC,"errmsg")
            end select
      end do

      do i = 1, nlocal_sondes
        n  = sonde_n(i)
        ix = sonde_x(i)
        iy = sonde_y(i)
!hf        n  = site_n(i)
!hf        ix = site_x(i)
!hf        iy = site_y(i)
        nn = 0


        !/** collect and print out with ground-level (KMAX_MID) first, hence &
        !    KMAX_MID:1:-1 in arrays  **/
        !/** first the advected and short-lived species

        do ispec = 1, NADV_SONDE    !/ xn_adv in ppb
           out(nn+1:nn+KMAX_MID,i) = PPBINV *  &
                                xn_adv( SONDE_ADV(ispec) , ix,iy,KMAX_MID:1:-1)
           nn = nn + KMAX_MID
        end do

        do ispec = 1, NSHL_SONDE    !/ xn_shl  in molecules/cm3 
           out(nn+1:nn+KMAX_MID,i) = xn_shl( SONDE_SHL(ispec) , ix,iy,KMAX_MID:1:-1)
           nn = nn + KMAX_MID
        end do

        !/** then print out XTRA stuff first, usually the height
        !    or pressure

        do ispec = 1, NXTRA_SONDE
          select case ( SONDE_XTRA(ispec) )
          case ( "z_mid" ) 
             out(nn+1:nn+KMAX_MID,i) =  z_mid(ix,iy,KMAX_MID:1:-1)   
          case ( "xksig" ) 
             out(nn+1:nn+KMAX_MID,i) =  xksig(ix,iy,KMAX_MID:1:-1)   ! NXTRA  first
          case ( "th" ) 
             out(nn+1:nn+KMAX_MID,i) =  th(ix,iy,KMAX_MID:1:-1,1)   ! NXTRA  first
          end select 
          nn = nn + KMAX_MID
        end do

      end do ! i


  !/** collect data into gout on me=0 t **/

   call siteswrt_out("sondes",IO_SONDES,NOUT_SONDE,NSONDES_MAX, FREQ_SONDE, &
                        nglobal_sondes,nlocal_sondes, &
                   sonde_gindex,sonde_name,sonde_gx,sonde_gy,sonde_species,out)

end subroutine siteswrt_sondes

!==================================================================== >

   subroutine siteswrt_out(fname,io_num,nout,nsites,f,nglobal,nlocal, &
                              s_gindex,s_name,s_gx,s_gy,s_species,out)

    !--- collect data from local nodes and writes out to sites/sondes.dat

    character(len=*), intent(in) :: fname
    integer, intent(in) :: io_num, nout, nsites
    integer, intent(in) :: f               ! Frequency of write-out (hours)
    integer, intent(in) :: nglobal, nlocal
    integer, intent(in), dimension (0:,:) :: s_gindex  ! index, starts at me=0
    character(len=*), intent(in), dimension (:) ::  s_name   ! site/sonde name
    integer, intent(in), dimension (:) :: s_gx, s_gy   ! coordinates
    character(len=*), intent(in), dimension (:) ::  s_species ! Variable names
    real,    intent(in), dimension(:,:) :: out    !/** outputs, local node **/

  ! Local
    real,dimension(nout,nglobal) :: g_out  !/** for output, collected  **/
    integer :: nglob, nloc, ix, iy               !/** Site indices           **/
    character(len=40) :: outfile
    character(len=4) :: suffix
    integer, parameter :: NTYPES = 2     ! No. types, now 2 (sites, sondes)
    integer ::  type                     ! = 1 for sites, 2 for sondes
    integer, save, dimension(NTYPES):: prev_month = (/ -99, -99 /) !Initialise

    select case (fname)
    case  ("sites" )
       type = 1
    case  ("sondes" )
       type = 2
    case default
       write(6,*) "non-possible tpye in siteswrt_out for ", fname
       return
    end select

    if(me == 0 .and. current_date%month /= prev_month(type) ) then

          if ( prev_month(type) > 0 ) close(io_num)  ! Close last-months file

         !/.. Open new file for write-out

          write(suffix,fmt="(2i2.2)") current_date%month, modulo ( current_date%year, 100 )
          outfile = fname // "." // suffix
          open(file=outfile,unit=io_num,action="write")
          prev_month(type) = current_date%month

          write(io_num,"(i3,2x,a,a, 4i4)") nglobal, fname, " in domain",  &
                                   gibegpos, giendpos, gjbegpos, gjendpos
          write(io_num,"(i3,a)") f, " Hours between outputs"

          do n = 1, nglobal
             write(io_num,"(a50,2i4)") s_name(n), s_gx(n), s_gy(n)
          end do ! nglobal

          write(io_num,"(i3,a)") size(s_species), " Variables:"
          do n = 1, size(s_species)
              write(io_num,"(i3,2x,a)") n, s_species(n)
          end do

   endif ! New month

  if ( me /= 0 ) then   ! -- send data to me=0

      call gc_isend(346,1,0,info,nloc,nlocal)
      if( nlocal > 0) call gc_rsend(347,nout*nlocal,0,info,out,out)

  else ! me = 0

      ! -- first, assign me=0 local data to g_out
      if ( MY_DEBUG ) print *, "ASSIGNS ME=0 NLOCAL_SITES", me, nlocal

      do n = 1, nlocal
         nglob = s_gindex(0,n)
         g_out(:,nglob) = out(:,n)
      end do ! n

      do d = 1, NPROC-1

          call gc_irecv(346,1,d,info,nloc,nloc)
          if( nloc > 0 ) call gc_rrecv(347,nout*nloc,d,info,out,out)

          do n = 1, nloc
             nglob = s_gindex(d,n)
             g_out(:,nglob) = out(:,n)
          end do ! n
      end do ! d


      ! Final output ***** 
      do n = 1, nglobal
          !MY_DEBUG write(6,*) ' Final out for n: ', n, s_name(n)
          write(unit=io_num, "(a10,i5,3i3,i5)" ) &
               s_name(n), current_date
          write(unit=io_num, "(5es10.3)" ) g_out(:,n)
      end do ! n

   end if ! me

  end subroutine siteswrt_out
end module Sites_ml
