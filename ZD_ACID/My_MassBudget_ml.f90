!_____________________________________________________________________________

module My_MassBudget_ml
!_____________________________________________________________________________
  use GenSpec_adv_ml      !! Can be many species
  use My_Emis_ml,     only :    QRCSO2,  QRCNO,  QRCNH3, QRCPM25, QRCPMco
  implicit none
  private

   !-----------------  "my" mass budget terms    ---------------------------!
   !  Here we define a few indices needed to relate species with IXADV_     !
   !  indices to their equivalent emission with QRC_ index.
   !                                                                        !  
   ! Plus, the array MY_MASS_PRINT to say which species to print out        !
   !------------------------------------------------------------------------!

   !-- contains subroutine:

      public :: set_mass_eqvs     ! Called from Emissions_ml

   ! Mass budget equivalency terms

    integer, public, parameter :: N_MASS_EQVS = 5!4  
    integer, public, save , dimension( N_MASS_EQVS ):: &
          ixadv_eqv  & !  IXADV_ no. of species
           ,qrc_eqv    !  QRC_   no. of equivalent species


   !/** species to print out in MassBudget_ml  (old myprint)
   !  Note - we can any number of species we need here - the dimensions 
   !  are obtained in MassBudget_ml with a size command.

   integer, public, parameter, dimension(12) :: MY_MASS_PRINT = &
     (/  IXADV_HNO3, IXADV_PAN, IXADV_NO, &
!hf amsu         IXADV_NO2,  IXADV_SO2, IXADV_SO4, IXADV_NH3, IXADV_AMSU, IXADV_AMNI/)
         IXADV_NO2,  IXADV_SO2, IXADV_SO4, IXADV_NH3, IXADV_aNH4, IXADV_aNO3 &
        ,IXADV_PM25, IXADV_PMco,  IXADV_pNO3  /)

  contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine set_mass_eqvs()
  !---------------------------------------------------------------------  
  !+  relates say IXADV_SO2 to QRCSO2

     ! Should have dimsnions N_MASS_EQVS

       ixadv_eqv(1) = IXADV_SO2
       ixadv_eqv(2) = IXADV_NO  
       ixadv_eqv(3) = IXADV_NH3 
       ixadv_eqv(4) = IXADV_PM25 
       ixadv_eqv(5) = IXADV_PMco

       qrc_eqv(1) = QRCSO2
       qrc_eqv(2) = QRCNO  
       qrc_eqv(3) = QRCNH3 
       qrc_eqv(4) = QRCPM25
       qrc_eqv(5) = QRCPMco


  end subroutine set_mass_eqvs
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module My_MassBudget_ml
