! <Runchem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                          module RunChem_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!----------------------------------------------------------------------
! Calls for routines calculating chemical and physical processes: 
! irreversible and equilibrium chemistry, dry and wet deposition,
! sea salt production, particle water etc.
!    
!----------------------------------------------------------------------
   use My_Aerosols_ml,    only: My_MARS, My_EQSAM, AERO_DYNAMICS,           &
                                EQUILIB_EMEP, EQUILIB_MARS, EQUILIB_EQSAM,  &
                                ORGANIC_AEROSOLS, Aero_water, SEASALT
   use My_Emis_ml       ! debug QRC
   use My_Timing_ml,      only: Code_timer, Add_2timing,  &
                                tim_before, tim_after

   use Ammonium_ml,       only: Ammonium
   use Aqueous_ml,        only: Setup_Clouds, prclouds_present, WetDeposition
   use Biogenics_ml      , only :  BIO_ISOP, BIO_TERP  ! debug only
   use CellMet_ml,        only: Get_CellMet
   use CheckStop_ml,      only: CheckStop
   use Chemfields_ml,     only: xn_adv  ! For DEBUG 
   use Chemsolver_ml,     only: chemistry
   use DefPhotolysis_ml,  only: setup_phot
   use DryDep_ml, only : drydep
   use ChemSpecs_tot_ml                   ! DEBUG ONLY
   use ChemSpecs_adv_ml                   ! DEBUG ONLY
   use GridValues_ml,     only: debug_proc, debug_li, debug_lj
   use ModelConstants_ml, only :  PPB, KMAX_MID, dt_advec, &
                                  nprint, END_OF_EMEPDAY, &
                  DebugCell, & ! DEBUG only
                  DEBUG => DEBUG_RUNCHEM, DEBUG_i, DEBUG_j,nstep, NPROC

   use OrganicAerosol_ml, only: OrganicAerosol ! not yet implemented 
   use Par_ml,            only : lj0,lj1,li0,li1  &
                                ,gi0, gj0, me & !! for testing
                                ,IRUNBEG, JRUNBEG    !! for testing
   use SeaSalt_ml,        only: SeaSalt_flux
   use Setup_1d_ml,       only: setup_1d, &
                                setup_bio, setup_rcemis, reset_3d
   use Setup_1dfields_ml, only: first_call  &
                     ,amk , rcemis, rcbio, xn_2d  ! DEBUG for testing
   use TimeDate_ml,       only: current_date

!--------------------------------
   implicit none
   private

   public :: runchem

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

subroutine runchem(numt)

   integer, intent(in) :: numt      

!/ local
   integer :: i, j
   integer :: errcode
   integer :: nmonth, nday, nhour     
   logical ::  Jan_1st, End_of_Run 
   logical ::  debug_flag    ! =>   Set true for selected i,j

! =============================

   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour

   Jan_1st    = ( nmonth == 1 .and. nday == 1 )
   End_of_Run = ( mod(numt,nprint) == 0       )

  !.... ****  processes calls *************************

   errcode = 0

    do j = lj0, lj1
      do i = li0, li1

          call Code_Timer(tim_before)

         !****** debug cell set here *******
          debug_flag =  .false.  
          debug_flag =  .false.  
          if ( DEBUG .and. debug_proc ) then

             debug_flag = ( debug_li == i .and. debug_lj == j ) 
             DebugCell = debug_flag
             if( debug_flag ) write(*,*) "RUNCHEM DEBUG START!"

          end if
         !write(*,"(a,4i4)") "RUNCHEM DEBUG IJTESTS", debug_li, debug_lj, i,j
         !write(*,*) "RUNCHEM DEBUG LLTESTS", me,debug_proc,debug_flag

 ! Prepare some near-surface grid and sub-scale meteorology
 ! for MicroMet
             call Get_CellMet(i,j,debug_flag) 

             call setup_1d(i,j)   

             call Add_2timing(27,tim_after,tim_before,&
                                              "Runchem:setup_1d")

             call Setup_Clouds(i,j,debug_flag)

             call setup_bio(i,j)    

             call Add_2timing(28,tim_after,tim_before,  &
                                         "Runchem:setup_cl/bio")

             call setup_phot(i,j,errcode)

             call CheckStop(errcode,"setup_photerror in Runchem") 
             call Add_2timing(29,tim_after,tim_before,  &
                                           "Runchem:1st setups")

             call setup_rcemis(i,j)

             if ( SEASALT )  &
             call SeaSalt_flux(i,j,debug_flag)

if ( DEBUG .and. debug_flag  ) then
    write(6,"(a,2i3,i5,9es10.2)") "DEBUG_RUNCHEM RCEMIS ", &
          current_date%day, current_date%hour, current_date%seconds, &
          rcemis(QRCNO,20), rcemis(QRCHCHO,20), rcemis(QRCISOP,20), rcbio(BIO_ISOP,20)
end if
!          !rcemis(QRCCO,20), AROM, rcemis(QRCPM25,20), rcemis(QRCEC_f_FFUEL,20)
!
!    write(6,"(a16,9es10.2)") "RUNCHEM PRE_OC ", &
!          xn_2d(PPM25,20),xn_2d(AER_BGNDOC,20), &
!          xn_2d(EC_F_FFUEL,20), xn_2d(POC_F_FFUEL,20), xn_2d(AER_FFUELOC,20),&
!          xn_2d(GAS_ASOA,20), xn_2d(AER_ASOA,20),&
!          xn_2d(GAS_BSOA,20), xn_2d(AER_BSOA,20)
!         
!          
             if ( ORGANIC_AEROSOLS )  &
               call OrganicAerosol(i,j,debug_flag)

             call Add_2timing(30,tim_after,tim_before,  &
                                          "Runchem:2nd setups")

!if ( DEBUG .and. debug_flag  ) then
!    write(6,"(a16,9es10.2)") "RUNCHEM POST_OC ", &
!          xn_2d(PPM25,20),xn_2d(AER_BGNDOC,20), &
!          xn_2d(EC_F_FFUEL,20), xn_2d(POC_F_FFUEL,20), xn_2d(AER_FFUELOC,20),&
!          xn_2d(GAS_ASOA,20), xn_2d(AER_ASOA,20),&
!          xn_2d(GAS_BSOA,20), xn_2d(AER_BSOA,20)
!end if
!
!
!           !-------------------------------------------------
!if ( DEBUG .and. debug_flag  ) then
!    write(6,"(a16,9es10.2)") "RUNCHEM PRE-CHEM ", &
!          xn_2d(PPM25,20), xn_2d(AER_BGNDOC,20), &
!          xn_2d(EC_F_FFUEL,20), xn_2d(POC_F_FFUEL,20), xn_2d(AER_FFUELOC,20),&
!          xn_2d(GAS_ASOA,20), xn_2d(AER_ASOA,20),&
!          xn_2d(GAS_BSOA,20), xn_2d(AER_BSOA,20)
!end if
             call chemistry(i,j)
!if ( DEBUG .and. debug_flag  ) then
!    write(6,"(a16,10es10.2)") "RUNCHEM POST-CHEM ", &
!          xn_2d(PPM25,20),xn_2d(AER_BGNDOC,20), &
!          xn_2d(EC_F_FFUEL,20), xn_2d(POC_F_FFUEL,20), xn_2d(AER_FFUELOC,20),&
!          xn_2d(GAS_ASOA,20), xn_2d(AER_ASOA,20),&
!          xn_2d(GAS_BSOA,20), xn_2d(AER_BSOA,20)
!end if
           !_________________________________________________

             call Add_2timing(31,tim_after,tim_before, &
                                           "Runchem:chemistry")
                
              !== Alternating Dry Deposition and Equilibrium chemistry
              !  Check that one and only one eq is chosen

                if(mod(nstep,2) /= 0 ) then 

                        if ( EQUILIB_EMEP )        call ammonium() 
                        if ( EQUILIB_MARS )        call My_MARS(debug_flag)
                        if ( EQUILIB_EQSAM )       call My_EQSAM(debug_flag) 

                        call DryDep(i,j)

                  else !do drydep first, then eq

                        call DryDep(i,j)
                        if ( EQUILIB_EMEP )        call ammonium() 
                        if ( EQUILIB_MARS )        call My_MARS(debug_flag)
                        if ( EQUILIB_EQSAM )       call My_EQSAM(debug_flag) 
                   endif
                   !????????????????????????????????????????????????????

                   call Add_2timing(32,tim_after,tim_before, &
                                                 "Runchem:ammonium+Drydep")

!                     if ( DEBUG .and. debug_flag  ) then
!                       write(6,"(a,10es12.3)") "DEBUG_RUNCHEM me RIEMER, aero", &
!                         xn_2d(SO4,20), xn_2d(aNO3,20), xn_2d(aNH4,20),tmpOut1, tmpOut2
!                       !write(6,*) "DEBUG_RUNCHEM me pre WetDep", me, prclouds_present
!                       write(6,"(a20,2i3,i5,3es12.3)") "DEBUG_RUNCHEM me OH", &
!                             current_date%day, current_date%hour,&
!                             current_date%seconds, &
!                             xn_2d(OH,20), xn_2d(O3,20), xn_2d(HNO3,20)
!
!                      end if

                     if ( prclouds_present)  &
                        call WetDeposition(i,j,debug_flag)

                   !** Modelling PM water at filter equlibration conditions:
                   !** T=20C and Rh=50% for comparability with gravimetric PM
 
                     if ( nhour == END_OF_EMEPDAY .or.  End_of_Run )    &
                        call Aero_water(i,j)                      

                     call reset_3d(i,j)

                     call Add_2timing(33,tim_after,tim_before,&
                                            "Runchem:post stuff")

                     first_call = .false.   !** end of first call **** !

      end do ! j
    end do ! i

                 !.... ************************************************

   end subroutine runchem

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module RunChem_ml

