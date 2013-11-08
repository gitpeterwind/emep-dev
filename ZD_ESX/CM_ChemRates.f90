!>_________________________________________________________<

module  ChemRates
!-----------------------------------------------------------

 
  use ChemFunctions       ! => kaero, RiemerN2O5

  use Zmet_ml        ! => tinv, h2o, m, Fgas !ESX
  use ZchemData      ! => rct
  use ChemSpecs      ! => NSPEC_TOT, PINALD, .... for FgasJ08
  use DefPhotolysis  ! => IDNO2 etc.
implicit none
  private

  !> CM_ChemRates RCTYPE rct sets Rate-coefficients - temperature dependant

    public :: setchemrates
    public :: setPhotolUsed

     integer, parameter, public :: NCHEMRATES = 70   !! No. coefficients

!> Photolysis rates
     integer, parameter, public :: NPHOTOLRATES = 13   !! No. DJ vals used
     integer, save, public,dimension(NPHOTOLRATES) :: photol_used
  contains
  !------------------------------------
  subroutine setPhotolUsed() 
    photol_used = (/ &
         IDO3_O3P &
        ,IDO3_O1D &
        ,IDNO2 &
        ,IDH2O2 &
        ,IDHNO3 &
        ,IDHCHO_H &
        ,IDHCHO_H2 &
        ,IDCH3CHO &
        ,IDNO3 &
        ,IDCH3O2H &
        ,IDCHOCHO &
        ,IDRCOCHO &
        ,IDMEK &
    /)

  end subroutine setPhotolUsed

  !------------------------------------
  subroutine setchemrates(debug_level) 
    integer, intent(in) :: debug_level

       rct(1,:) = (6.0e-34*O2+5.6E-34*N2)*O2*exp(-2.6*LOGTDIV300) 
       rct(2,:) = 1.8e-11*N2*exp(107.0*TINV) 
       rct(3,:) = 3.2e-11*O2*exp(67.0*TINV) 
       rct(4,:) = 2.2e-10*H2O 
       rct(5,:) = 1.4e-12*exp(-1310.0*TINV) 
       rct(6,:) = 1.4e-13*exp(-2470.0*TINV) 
       rct(7,:) = 1.7e-12*exp(-940.0*TINV) 
       rct(8,:) = 2.03e-16*exp(-4.57*LOG300DIVT)*exp(693.0*TINV) 
       rct(9,:) = 1.8e-11*exp(110.0*TINV) 
       rct(10,:) = 3.6e-12*exp(270.0*TINV) 
       rct(11,:) = 4.5e-14*exp(-1260.0*TINV) 
       rct(12,:) = 4.8e-11*exp(250.0*TINV) 
       rct(13,:) = 2.9e-12*exp(-160.0*TINV) 
       rct(14,:) = 7.7e-12*exp(-2100.0*TINV) 
       rct(15,:) = KMT3(2.4e-14,460.0,6.5E-34,1335.0,2.7E-17,2199.0) 
       rct(16,:) = (1.0+1.4e-21*H2O*exp(2200.0*TINV))*2.2E-13*exp(600.0*TINV) 
       rct(17,:) = (1.0+1.4e-21*H2O*exp(2200.0*TINV))*1.9E-33*exp(980.0*TINV)*M 
       rct(18,:) = 2.5e-12*exp(-260.0*TINV) 
       rct(19,:) = 1.85e-20*exp(2.82*LOG(TEMP))*exp(-987.0*TINV) 
       rct(20,:) = 1.44e-13+M*3.43E-33 
       rct(21,:) = 2.3e-12*exp(360.0*TINV) 
       rct(22,:) = 7.4e-13*exp(-520.0*TINV) 
       rct(23,:) = 1.03e-13*exp(365.0*TINV)-7.4E-13*exp(-520.0*TINV) 
       rct(24,:) = 6.38e-18*(TEMP**2)*exp(144.0*TINV) 
       rct(25,:) = 3.8e-13*exp(780.0*TINV) 
       rct(26,:) = 5.3e-12*exp(190.0*TINV) 
       rct(27,:) = 1.25e-17*(TEMP**2)*exp(615.0*TINV) 
       rct(28,:) = 2e-12*exp(-2440.0*TINV) 
       rct(29,:) = 6.9e-12*exp(-1000.0*TINV) 
       rct(30,:) = 2.55e-12*exp(380.0*TINV) 
       rct(31,:) = 3.8e-13*exp(900.0*TINV) 
       rct(32,:) = 1.9e-12*exp(190.0*TINV) 
       rct(33,:) = 4.4e-12*exp(365.0*TINV) 
       rct(34,:) = 7.5e-12*exp(290.0*TINV) 
       rct(35,:) = 2e-12*exp(500.0*TINV) 
       rct(36,:) = 2.9e-12*exp(500.0*TINV) 
       rct(37,:) = 5.2e-13*exp(980.0*TINV) 
       rct(38,:) = 6.7e-18*(TEMP**2)*exp(511.0*TINV) 
       rct(39,:) = 2.03e-17*(TEMP**2)*exp(78.0*TINV) 
       rct(40,:) = 2.54e-12*exp(360.0*TINV) 
       rct(41,:) = 1.81875e-13*exp(1300.0*TINV) 
       rct(42,:) = 2.53e-18*(TEMP**2)*exp(503.0*TINV) 
       rct(43,:) = 9.1e-15*exp(-2580.0*TINV) 
       rct(44,:) = 5.5e-15*exp(-1880.0*TINV) 
       rct(45,:) = 1.5132e-13*exp(1300.0*TINV) 
       rct(46,:) = 2.49969e-13*exp(1300.0*TINV) 
       rct(47,:) = 2.05446e-13*exp(1300.0*TINV) 
       rct(48,:) = 6.6e-18*(TEMP**2)*exp(820.0*TINV) 
       rct(49,:) = 1.9e-12*exp(575.0*TINV) 
       rct(50,:) = 1.03e-14*exp(-1995.0*TINV) 
       rct(51,:) = 2.7e-11*exp(390.0*TINV) 
       rct(52,:) = 2.6e-12*exp(610.0*TINV) 
       rct(53,:) = 1.36e-15*exp(-2112.0*TINV) 
       rct(54,:) = 8e-12*exp(380.0*TINV) 
       rct(55,:) = 7.6e-12*exp(180.0*TINV) 
       rct(56,:) = 1.6e-12*exp(305.0*TINV) 
       rct(57,:) = 8.7e-12*exp(290.0*TINV) 
       rct(58,:) = 8.5e-16*exp(-1520.0*TINV) 
       rct(59,:) = 3.15e-12*exp(-450.0*TINV) 
       rct(60,:) = 4.3e-13*exp(1040.0*TINV) 
       rct(61,:) = KAERO() 
       rct(62,:) = IUPAC_TROE(1.0e-31*exp(1.6*LOG300DIVT)  &
         ,3.0E-11*exp(-0.3*LOG300DIVT)  &
         ,0.85  &
         ,M  &
         ,0.75-1.27*LOG10(0.85)) 
       rct(63,:) = IUPAC_TROE(3.6e-30*exp(4.1*LOG300DIVT)  &
         ,1.9E-12*exp(-0.2*LOG300DIVT)  &
         ,0.35  &
         ,M  &
         ,0.75-1.27*LOG10(0.35)) 
       rct(64,:) = IUPAC_TROE(1.3e-3*exp(3.5*LOG300DIVT)*exp(-11000.0*TINV)  &
         ,9.70E14*exp(-0.1*LOG300DIVT)*exp(-11080.0*TINV)  &
         ,0.35  &
         ,M  &
         ,0.75-1.27*LOG10(0.35)) 
       rct(65,:) = IUPAC_TROE(3.3e-30*exp(3.0*LOG300DIVT)  &
         ,4.1E-11  &
         ,0.40  &
         ,M  &
         ,0.75-1.27*LOG10(0.4)) 
       rct(66,:) = IUPAC_TROE(2.7e-28*exp(7.1*LOG300DIVT)  &
         ,1.2E-11*exp(0.9*LOG300DIVT)  &
         ,0.3  &
         ,M  &
         ,0.75-1.27*LOG10(0.3)) 
       rct(67,:) = IUPAC_TROE(4.9e-3*exp(-12100.0*TINV)  &
         ,5.4E16*exp(-13830.0*TINV)  &
         ,0.3  &
         ,M  &
         ,0.75-1.27*LOG10(0.3)) 
       rct(68,:) = IUPAC_TROE(8.6e-29*exp(3.1*LOG300DIVT)  &
         ,9.0E-12*exp(0.85*LOG300DIVT)  &
         ,0.48  &
         ,M  &
         ,0.75-1.27*LOG10(0.48)) 
       rct(69,:) = IUPAC_TROE(8.0e-27*exp(3.5*LOG300DIVT)  &
         ,3.0E-11*300.0*TINV  &
         ,0.5  &
         ,M  &
         ,0.75-1.27*LOG10(0.5)) 
       rct(70,:) = IUPAC_TROE(7.4e-31*exp(2.4*LOG300DIVT)  &
         ,3.3E-11*exp(0.3*LOG300DIVT)  &
         ,exp(-temp/1420.0)  &
         ,M  &
         ,0.75+3.884E-4*temp) 

  end subroutine setchemrates
end module  ChemRates
