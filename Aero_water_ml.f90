      module Aero_water_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   calculates aerosols liquid water content
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 implicit none
 private
!-- subroutines:
        public ::  Awater

!-- functions
	private ::  poly4, poly6

   contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine Awater(relh,mso4,mnh4,mno3,wh2o)
  !.......................................................
! note!!! wh2o is returned in micrograms / cubic meter
!         mso4,mnh4,mno3 are in micromoles / cubic meter
!
!  this  version uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of water
! activity where:
!                 mfs = ms / ( ms + mw)
!                 ms - the mass of solute
!                 mw - the mass of water.
!
!  define   y = mw/ ms
!
!  then     mfs = 1 / (1 + y)
!
!    y can then be obtained from the values of mfs as
!
!             y = (1 - mfs) / mfs
!
!     The aerosol is assumed to be in a metastable state if the rh is
!     is below the rh of deliquescence, but above the rh of crystallization.
!
!     zsr interpolation is used for sulfates with x ( the molar ratio of
!     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
!     section 1: 0 <= x < 1
!     section 2: 1 <= x < 1.5
!     section 3: 1.5 <= x < 2.0
!     section 4: 2 <= x
!     In sections 1 through 3, only the sulfates can affect the amount of water
!     on the particles.
!     In section 4, we have fully neutralized sulfate, and extra ammonium which
!     allows more nitrate to be present. thus, the ammount of water is calculated
!     using zsr for ammonium sulfate and ammonium nitrate. Crystallization is
!     assumed to occur in sections 2,3,and 4. see detailed discussion below.
!
! Definitions:
!     mso4, mnh4, and mno3 - the number of micromoles/(cubic meter of air)
!      for sulfate, ammonium, and nitrate respectively,
!     irhx - the relative humidity (%)
!     wh2o - the returned water amount in micrograms / cubic meter of air
!     x - the molar ratio of ammonium to sulfate
!     y0,y1,y1.5, y2 - the water contents in mass of water/mass of solute
!     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
!     y3 - the value of the mass ratio of water to solute for
!     a pure ammonium nitrate  solution.
!
! .. Coded by dr. francis s. binkowski, 4/8/96.
!.........................................................................

 use GenSpec_tot_ml,  only : SO4, aNH4, aNO3        !  For mol. wts.
 use GenChemicals_ml, only : species               !  For mol. wts.

!hf use GenSpec_aer_ml,  only : MAITSO4, MAITNO3, MAITNH4     !  For mol. wts.
!hf use GenChem_aer_ml,  only : aerospec

implicit none

  !-- in          
    real, intent(in)    :: relh, mso4, mnh4, mno3
  !-- out
    real, intent(out)  ::  wh2o

  !-- local
    real   ::  tso4,tnh4,tno3, x, awc, u, irh , aw              &
              !,poly4, poly6, mfs0,mfs1,mfs15, mfs2         &
              , mfs0,mfs1,mfs15, mfs2         &
              ,y, y0,y1,y15,y2,y3,y40,y140,y1540,yc        &
              , mfsso4, mfsno3
    real, dimension(4) :: C0, C1, C15, C2           
    real, dimension(6) ::  KSO4, KNO3
      real   ::   mwso4,mwnh4,mwno3, mw2, mwano3      

 !   real, parameter ::                   &
            !  ,  mwso4   = 96.0636      &
            !  ,  mwnh4   = 18.0985      &
            !  ,  mwno3   = 62.0649      
 
!-------------------------------------------------------------------
!     the polynomials use data for relh as a function of mfs from Tang and
!     Munkelwitz, JGR. 99: 18801-18808, 1994.
!     The polynomials were fit to tang's values of water activity as a
!     function of mfs.

! *** coefficients of polynomials fit to Tang and Munkelwitz data
!     now give mfs as a function of water activity.
!---------------------------------------------------------------
      data C1/0.9995178, -0.7952896, 0.99683673, -1.143874/
      data C15/1.697092,-4.045936, 5.833688, -3.463783/
      data C2/2.085067, -6.024139, 8.967967, -5.002934/
!---------------------------------------------------------------
! *** the following coefficients are a fit to the data in table 1 of
!     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
!      data C0/0.8258941, -1.899205, 3.296905, -2.214749 /
! *** new data fit to data from
!       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
!       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
!       Zeleznik J. Phys. Chem. ref. data, 20: 157-1200
!---------------------------------------------------------------
       data C0/ 0.798079, -1.574367, 2.536686, -1.735297 /
!---------------------------------------------------------------
! *** polynomials for ammonium nitrate and ammonium sulfate are from:
!     chan et al.1992, atmospheric environment (26a): 1661-1673.
!---------------------------------------------------------------
      data KNO3/0.2906, 6.83665, -26.9093, &
     &          46.6983, -38.803, 11.8837/
      data KSO4/ 2.27515, -11.147, 36.3369, &
     &        -64.2134, 56.8341, -20.0953/
!---------------------------------------------------------------

!//.. check range of per cent relative humidity
!     aw - water activity = fractional relative humidity

       aw = relh
       aw = max(0.01,aw)
! ia actionia TEST SET MAXIMUM TO 97
       aw = min(aw,0.95)	! modified for consistency 100->99 IA

!====================================================================

!hf         mwso4   = aerospec(MAITSO4)%molwt
!hf         mwnh4   = aerospec(MAITNH4)%molwt      
!hf         mwno3   = aerospec(MAITNO3)%molwt

            mwso4   = species(SO4)%molwt
            mwnh4   = species(aNH4)%molwt
            mwno3   = species(aNO3)%molwt

         mw2  = mwso4 + 2.0 * mwnh4   
         mwano3  = mwno3 + mwnh4

       tso4 = max( mso4 , 0.0 )
       tnh4 = max( mnh4 , 0.0 )
       tno3 = max( mno3 , 0.0 )

       x = 0.0
! *** if there is non-zero sulfate calculate the molar ratio
       if (tso4 > 0.0 ) then
             x = tnh4 / tso4
       else
! *** otherwise check for non-zero nitrate and ammonium
         if ( tno3 > 0.0 .and. tnh4 > 0.0 )   x = 10.0
       end if
!
! *** begin screen on x for calculating wh2o
       if ( x < 1.0 ) then
!
          mfs0 = poly4(C0,aw)
          mfs1 = poly4(C1,aw)
          y0 = (1.0 - mfs0 ) / mfs0
          y1 = (1.0 - mfs1 ) / mfs1
          y = (1.0 - x) * y0 + x * y1

!
       else if ( x < 1.5) then
!
         if ( aw >= 0.40 ) then
            mfs1  = poly4(C1,aw)
            mfs15 = poly4(C15,aw)
            y1  = (1.0 - mfs1 ) / mfs1
            y15 = (1.0 - mfs15) / mfs15
            y = 2.0 * ( y1 * (1.5 - x) + y15 *( x - 1.0) )
         else
! *** set up for crystalization

! *** crystallization is done as follows:
!      for 1.5 <= x, crystallization is assumed to occur at rh = 0.4
!      for x <= 1.0, crystallization is assumed to occur at an rh < 0.01,
!      and since the code does not allow ar rh < 0.01, crystallization
!      is assumed not to occur in this range.
!      for 1.0 <= x <= 1.5 the crystallization curve is a straignt line
!      from a value of y15 at rh = 0.4 to a value of zero at y1. from
!      point b to point a in the diagram.
!      the algorithm does a double interpolation to calculate the amount of
!      water.
!
!        y1(0.40)               y15(0.40)
!         +                     + point b
!
!
!
!
!         +--------------------+
!       x=1                   x=1.5
!      point a
!
!

            awc = 0.80 * (x - 1.0) ! rh along the crystallization curve.

            y = 0.0
               u=0.40
             if ( aw >= awc ) then    ! interpolate using crystalization curve
               mfs1  = poly4(C1,u)
               mfs15 = poly4(C15,u)
               y140  = (1.0 - mfs1 ) / mfs1
               y1540 = (1.0 - mfs15) / mfs15
               y40 = 2.0 * ( y140 * (1.5 - x) + y1540 *( x - 1.0) )
               yc = 2.0 * y1540 * (x -1.0) ! y along crystallization curve
               y = y40 - (y40 - yc) * (u - aw) / (u - awc)
            end if ! end of checking for aw
          end if ! end of checking on irh

       else if( x < 1.9999) then
!
           y= 0.0
           if( aw >= 0.40) then
             mfs15 = poly4(C15,aw)
             mfs2  = poly4(C2,aw)
             y15 = (1.0 - mfs15) / mfs15
             y2  = (1.0 - mfs2) / mfs2
             y = 2.0 * (y15 * (2.0 - x) + y2 * (x - 1.5) )
           end if ! end of check for crystallization
!
      else          !  1.9999 < x

! regime where ammonium sulfate and ammonium nitrate are in solution.
!
! *** following cf&s for both ammonium sulfate and ammonium nitrate
! *** check for crystallization here. their data indicate a 40% value
!     is appropriate.
            y2 = 0.0
            y3 = 0.0
            if ( aw >= 0.40) then
              mfsso4 = poly6(KSO4,aw)
              mfsno3 = poly6(KNO3,aw)
              y2 = (1.0 - mfsso4) / mfsso4
              y3 = (1.0 - mfsno3) / mfsno3

            end if
!
       end if ! end of checking on x
!
! *** now set up output of wh2o

!      wh2o units are micrograms (liquid water) / cubic meter of air
!
       if ( x < 1.9999) then

         wh2o =  y * (tso4 * mwso4 + tnh4 * mwnh4 )

       else

! *** this is the case that all the sulfate is ammonium sulfate
!     and the excess ammonium forms ammonum nitrate

        wh2o =   y2 * tso4 * mw2 + y3 * tno3 * mwano3

       end if


       end subroutine Awater
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      real function poly4 (a,x)
  !.......................................
    implicit none

  !-- arguments
    real, dimension(4), intent(in)  ::  a
    real, intent(in)  ::  x
!------------------------------------------------------------

       poly4 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) )))
      
      return 
      end function poly4
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      real function poly6(a,x)
  !.......................................
    implicit none

  !-- arguments
    real, dimension(6), intent(in)  ::  a
    real, intent(in)  ::  x
!------------------------------------------------------------
      poly6 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) +     &
     &           x * ( a(5) + x * (a(6)  )))))

      return       
      end  function poly6 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   

      end module Aero_water_ml
