!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 module   MassBudget_ml

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! DESCRIPTION 
! ........ to be added !!!
! 1/10/01 - code for derived fields removed. MY_MASS_PRINT ADDED, ds
! Oct, 2001 - ** new ** mass budget method by jej
! Nov.2001, ds, "use" statements moved to top, much code moved to REMOVED 
! section at end
!_____________________________________________________________________________

!hf u2 use My_Runmode_ml , only : DEBUG
!ds rv1.2 use My_Derived_ml,  only : NWDEP, NDDEP & ! No. deposition fields
!ds rv1.2                        ,f_wdep, f_ddep  & ! definitions of dep data fields
!ds rv1.2                        , wdep, ddep       ! deposition data fields
 use My_MassBudget_ml,only: MY_MASS_PRINT  ! Species to be printed
                                           ! (old myprint array)

 use GenChemicals_ml, only: species
 use GenSpec_adv_ml, only : NSPEC_ADV      ! No. species (long-lived)
 use GenSpec_shl_ml, only : NSPEC_SHL      ! No. species (shorshort-lived)
 !u1 use GenSpec_maps_ml, only: MAP_ADV2TOT     ! Index mapping
 use Chemfields_ml , only : xn_adv
 use GridValues_ml , only : carea,xmd
 use Io_ml         , only : IO_RES
 use Met_ml        , only : ps, psurf   !u7.4vg - was psa
 use ModelConstants_ml , &
                     only : KMAX_MID  &  ! Number of points (levels) in vertical
                           ,PT        &  ! Pressure at top
                           ,ATWAIR, atwS, atwN
 use Par_ml,only:  & 
        MAXLIMAX,  & ! Maximum number of local points in longitude
        MAXLJMAX,  & ! Maximum number of local points in latitude
        li0,li1,lj0,lj1,NPROC,me, limax, ljmax, &
        gi0, gj0, GIMAX, GJMAX
 
implicit none
private


! The following parameters are used to check the global mass budget:
! Initialise here also.

  real, public, save, dimension(NSPEC_ADV) ::   &
      sumint   = 0.0   & !  initial mass
     ,fluxin   = 0.0   & !  mass in  across lateral boundaries
     ,fluxout  = 0.0   & !  mass out across lateral boundaries
     ,totddep  = 0.0   & !  total dry dep
     ,totwdep  = 0.0   & !  total wet dep
     ,totem    = 0.0   & !  total emissions
     ,totox    = 0.0   & !  total oxidation
     ,totldep  = 0.0     !  local deposition (not in use - Lagrangian)

  real, public, save, dimension(NSPEC_ADV) ::  &
      amax = -2.0   &  ! maximum concentration in field -2
     ,amin =  2.0      ! minimum concentration in field  2

!hf u2:
  logical, private, parameter :: DEBUG = .false.

  public :: Init_massbudget
  public :: massbudget
  

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine Init_massbudget()
    ! Initialise mass-budget - calculate mass of concentrations fields
    ! within 3-D grid, after boundary conditions
    !
    ! COnverted from old tsfld , ds, 14/5/01
    !-------------------------------------------------------------------------

    integer i, j, k, n, info
    real rwork

    do k=2,KMAX_MID   
      do j=lj0,lj1
        do i=li0,li1
            rwork = carea(k)* xmd(i,j)*(ps(i,j,1) - PT)
            sumint(:) = sumint(:) + xn_adv(:,i,j,k)*rwork  ! sumint in kg      
        enddo
      enddo
    enddo

    call gc_rsum(NSPEC_ADV, NPROC, info, sumint)
    if(me == 0)then
         do n = 1,NSPEC_ADV
	   if(sumint(n) >  0. ) then
             write(IO_RES,"(a15,i2,4x,e10.3)") "Initial mass",n,sumint(n) 
             write(6,"(a15,i2,4x,e10.3)") "Initial mass",n,sumint(n) 
           end if
         enddo
    end if

 end subroutine Init_massbudget


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine massbudget()

  ! Converted from old masbud.f by ds, March 2001
  !   sums over all sulphur and nitrogen, so is model independant.


  integer ::  i, j, k, n, nn, info    !6b, lpr1, lpr2
  integer ispec
  real, dimension(NSPEC_ADV,KMAX_MID) ::  sumk   ! total mass in each layer
  character(len=12) :: spec_name
  integer ispec_name

!pw  real, dimension(NSPEC_ADV) ::  &
!pw           amax, amin        ! max/min of conc. (long lived)

  real, dimension(3)         :: family_int  & ! initial total mass of 
                                              ! species family
                               ,family_mass & ! total family mass at the 
                                              ! end of the model run
                               ,family_in  &  ! total family mass flowing in  
                               ,family_out &  ! total family mass flowing out 
                               ,family_ddep&  ! total family mass dry dep. 
                               ,family_wdep&  ! total family mass wet dep. 
                               ,family_em  &  ! total family mass emitted 
                               ,family_div &  ! total family mass input
                               ,family_frac   ! total family mass output 

  real, dimension(NSPEC_ADV) :: xmax, xmin, & ! min and max value for the 
					      ! individual species
                                sum_mass,   & ! total mass of species
                                frac_mass,  & ! mass budget fraction (should 
					      ! be one) for groups of species
     			gfluxin, gfluxout,  & ! flux in  and out
				 gtotem,    & ! total emission
			gtotddep, gtotwdep, & ! total dry and wet deposition
				  gtotldep, & ! local dry deposition
				  gtotox      ! oxidation of SO2  ????????


  real :: totdiv,helsum,ammfac, atw, natoms

    sum_mass(:)     = 0.
    frac_mass(:)    = 0.
    xmax(:)         = -2.
    xmin (:)        = 2.
    gfluxin(:)   = fluxin(:)
    gfluxout(:)  = fluxout(:)
    gtotem(:)    = totem(:)
    gtotddep(:)  = totddep(:)
    gtotwdep(:)  = totwdep(:)
    gtotldep(:)  = totldep(:)
    gtotox(:)    = totox(:)

    sumk(:,:) = 0.

    do k = 1,KMAX_MID
      do j = lj0,lj1
        do i = li0,li1
!            helsum = carea(k)*xmd(i,j) * (ps(i,j,1) - PT)
            helsum = carea(k)*xmd(i,j) * (psurf(i,j) - PT)

            xmax(:) = amax1(xmax(:),xn_adv(:,i,j,k))
            xmin(:) = amin1(xmin(:),xn_adv(:,i,j,k))

            sumk(:,k) = sumk(:,k) + xn_adv(:,i,j,k)*helsum

        enddo
      enddo
    enddo



    call gc_rmax(NSPEC_ADV, NPROC, info, xmax)
    call gc_rmin(NSPEC_ADV, NPROC, info, xmin)
    call gc_rsum(NSPEC_ADV, NPROC, info, gfluxin)
    call gc_rsum(NSPEC_ADV, NPROC, info, gfluxout)
    call gc_rsum(NSPEC_ADV, NPROC, info, gtotem)
    call gc_rsum(NSPEC_ADV, NPROC, info, gtotddep)
    call gc_rsum(NSPEC_ADV, NPROC, info, gtotwdep)
    call gc_rsum(NSPEC_ADV, NPROC, info, gtotldep)
    call gc_rsum(NSPEC_ADV, NPROC, info, gtotox)
    call gc_rsum(NSPEC_ADV*KMAX_MID, NPROC, info, sumk)

!     make some temporary variables used to hold the sum over all 
!     domains. Remember that sumint already holds the sum over all
!     domains, see inass
!
    amax(:) = max( amax(:), xmax(:) )
    amin(:) = min( amin(:), xmin(:) )
    do k = 2,KMAX_MID
      sum_mass(:) = sum_mass(:)+sumk(:,k)
    enddo



    do n = 1,NSPEC_ADV
!      totdiv = sumint(n) + gtotem(n)* &
!               ATWAIR/species( MAP_ADV2TOT(n))%molwt + gfluxin(n)
      totdiv = sumint(n) + gtotem(n) + gfluxin(n)
      frac_mass(n) = sum_mass(n)  + (gtotddep(n)+gtotwdep(n))*ATWAIR &
                   + gfluxout(n) 

!    NO LOCAL DEPOSITION DEFINED
!                   + (gtotldep(n))*ATWAIR/species( MAP_ADV2TOT(n))%molwt &
!!  
      if(totdiv >  0.0 ) frac_mass(n) = frac_mass(n)/totdiv


    end do


   if (me == 0 ) then

    do n = 1,NSPEC_ADV
      if (gtotem(n) > 0.0 ) write(6,*)          &
                           'tot. emission of species ',n,gtotem(n)
    end do

    family_int(:)  = 0.
    family_mass(:) = 0.
    family_in (:)  = 0.
    family_out(:)  = 0.
    family_div(:)  = 0.
    family_frac(:) = 0.
    family_ddep(:) = 0.
    family_wdep(:) = 0.
    family_em(:)   = 0.

     do ispec = 1, 3

       write(6,*) 'family ', ispec 
       do n = 1, NSPEC_ADV

         !u1 nn = MAP_ADV2TOT(n)
         nn = NSPEC_SHL + n


           if ( ispec == 1 ) natoms = real(species(nn)%sulphurs)
           if ( ispec == 1 ) atw = atwS

           if ( ispec == 2 ) natoms = real(species(nn)%nitrogens)
           if ( ispec == 2 ) atw = atwN

           if ( ispec == 3 ) natoms = real(species(nn)%carbons)
           if ( ispec == 3 ) atw = 16
         
           if (natoms > 0) then
             family_int(ispec) = family_int(ispec) + sumint(n)*natoms
             family_mass(ispec) = family_mass(ispec) + sum_mass(n)*natoms
             family_in (ispec) = family_in (ispec) + gfluxin(n)*natoms
             family_out(ispec) = family_out(ispec) + gfluxout(n)*natoms
             family_ddep(ispec) = family_ddep(ispec) + gtotddep(n)*natoms
             family_wdep(ispec) = family_wdep(ispec) + gtotwdep(n)*natoms
             family_em(ispec) = family_em(ispec) + gtotem(n)*natoms
           end if
       end do  ! NSPEC_ADV

       family_div(ispec) = family_int(ispec) &
                         + family_in (ispec) &
                         + family_em(ispec)   !6b &
!!                        * ATWAIR/atw
!                         * ATWAIR/species( MAP_ADV2TOT(n))%molwt

       if (family_div(ispec) > 0.0 ) &
              family_frac(ispec) = (family_mass(ispec) &
                                 +  family_out(ispec)  &
                                 +  family_ddep(ispec)*ATWAIR  & ! not /atw(n)
                                 +  family_wdep(ispec)*ATWAIR) & ! not /atw(n)
                                 / family_div(ispec)


      write(6,*)
      write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'      
      write(6,*)

      !u4 write(6,951)
      !u4 write(6,952) ispec,family_int(ispec), family_mass(ispec) &

      write(6,"(6a12)") "parameter", "sumint", "summas", &
                               "fluxout","fluxin", "fracmass"
      write(6,"(i9,5es12.4)") ispec,family_int(ispec), family_mass(ispec) &
                        ,family_out(ispec), family_in (ispec)  &
                        ,family_frac(ispec)

      write(6,*)
      write(6,958) 
      write(6,957) ispec, family_ddep(ispec)*ATWAIR  &
                        , family_wdep(ispec)*ATWAIR  &
                        , family_em(ispec) 
      write(6,*)
      write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'

    end do  ! ispec = 1,3
    

   end if

  if (me == 0) then     ! printout from node 0

     !/.. now use species array which is set in My_MassBudget_ml
    do k = 1,KMAX_MID

      write(6,*)
      write(IO_RES,*)

      do nn = 1,size ( MY_MASS_PRINT )
        n = MY_MASS_PRINT(nn)
        write(6,950)      n,k,sumk(n,k)
        write(IO_RES,950) n,k,sumk(n,k)
      end do
    enddo


     do nn = 1,size ( MY_MASS_PRINT )
        n = MY_MASS_PRINT(nn)

        write(6,*)
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'      
        write(6,*)
        !u4 write(6,951)
        !u4 write(6,952) n,sumint(n),sum_mass(n),gfluxout(n),gfluxin(n), &

        write(6,"(6a12)") "parameter", "sumint", "summas", &
                               "fluxout","fluxin", "fracmass"
        write(6,"(i9,5es12.4)") n,sumint(n),sum_mass(n), &
                      gfluxout(n),gfluxin(n), frac_mass(n)

        write(6,*)
        write(6,955) 
        write(6,956) n, gtotox(n), gtotddep(n), gtotwdep(n),         &
                   gtotem(n), gtotldep(n)
        write(6,*)
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
                 
    enddo
!              





   end if  ! me = 0

 950  format(1h,'parameter ',i2,5x,'level',i2,5x,es12.5)
!u4 951  format(1h,'parameter',7x,'sumint',8x,'summas',8x,'fluxout ',8x, &
!u4                'fluxin  ',6x,'fracmas')
!u4 952  format(1h ,6x,i2,7x,es12.5,4x,es12.5,4x,es12.5,4x,es12.5,4x,es12.5)
 953  format(1h ,9x,'amax  ',8x,'amin  ',8x,'xmax  ',8x,'xmin  ')
 954  format(1h ,i2,4x,e10.3,4x,e10.3,4x,e10.3,4x,e10.3)
 955  format(1h ,9x,'totox',16x,'totddep',13x,                        &
          'totwdep',13x,'totem',13x,'totldep')
 956  format(1h ,i2,4x,es10.3,10x,es10.3,10x,es10.3,10x,es10.3,10x,es10.3)
 957  format(1h ,i2,4x,es10.3,10x,es10.3,10x,es10.3)
 958  format(1h ,9x,'totddep',13x,'totwdep',13x,'totem')

 end subroutine massbudget


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
             end module MassBudget_ml
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
