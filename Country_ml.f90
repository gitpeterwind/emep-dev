module Country_ml

 !+ Sets country index numbers (IC_xx), code, time-zones, and names
 !
 ! Language : F
 ! History  :
 !   xnd version: 13/11/2003 - jej/ds, new countries, corrected 29
 !   2nd version: 12/12/2001 - ds,  added is_sea to country type
 !   1st version: 24/1/2001 - ds, using timezones from Vigdis.
 !
 ! notes: Unfortunately doing the setting of the Country data needs the
 ! subroutine Set_Countries, as doing this with a neat "parameter"
 ! statement didn't work - there were too many continuation lines :-(
 ! And I would have used the word countries instead of cc, but that made
 ! the lines too long in Set-Countries....

  ! In this module we give details for all countries in the UNECE tables,
  ! except that the US, Canada and EU are commented out.
  implicit none

  public :: Country_Init     ! sets country details

  integer, parameter, public :: NLAND = 71
  logical, parameter, private :: T = .true.   ! shorthand
  logical, parameter, private :: F = .false.  ! shorthand

  !/ to be set in Country_Init:

  type, public :: cc
     character(len=3)  :: code
     integer           :: index
     logical           :: is_sea     ! 1 for sea area, 0 otherwise
     integer           :: timezone
     character(len=30) :: name
  end type cc

  type(cc), public, save, dimension(NLAND) :: Country

  integer, parameter, public ::  IC_AL =   1   ! Albania                       
  integer, parameter, public ::  IC_AT =   2   ! Austria                       
  integer, parameter, public ::  IC_BE =   3   ! Belgium                       
  integer, parameter, public ::  IC_BG =   4   ! Bulgaria                      
  integer, parameter, public ::  IC_CS =   5   ! Former Yugoslavia             
  integer, parameter, public ::  IC_DK =   6   ! Denmark                       
  integer, parameter, public ::  IC_FI =   7   ! Finland                       
  integer, parameter, public ::  IC_FR =   8   ! France                        
  integer, parameter, public :: IC_GDR =   9   ! Former East Germany           
  integer, parameter, public :: IC_FRG =  10   ! Former West Germany           
  integer, parameter, public ::  IC_GR =  11   ! Greece                        
  integer, parameter, public ::  IC_HU =  12   ! Hungary                       
  integer, parameter, public ::  IC_IS =  13   ! Iceland                       
  integer, parameter, public ::  IC_IE =  14   ! Ireland                       
  integer, parameter, public ::  IC_IT =  15   ! Italy                         
  integer, parameter, public ::  IC_LU =  16   ! Luxembourg                    
  integer, parameter, public ::  IC_NL =  17   ! Netherlands                   
  integer, parameter, public ::  IC_NO =  18   ! Norway                        
  integer, parameter, public ::  IC_PL =  19   ! Poland                        
  integer, parameter, public ::  IC_PT =  20   ! Portugal                      
  integer, parameter, public ::  IC_RO =  21   ! Romania                       
  integer, parameter, public ::  IC_ES =  22   ! Spain                         
  integer, parameter, public ::  IC_SE =  23   ! Sweden                        
  integer, parameter, public ::  IC_CH =  24   ! Switzerland                   
  integer, parameter, public ::  IC_TR =  25   ! Turkey                        
  integer, parameter, public ::  IC_SU =  26   ! Former USSE                   
  integer, parameter, public ::  IC_GB =  27   ! United Kingdom                
  integer, parameter, public :: IC_VUL =  28   ! Vulcanoes                     
  integer, parameter, public :: IC_REM =  29   ! Remaining                     
  integer, parameter, public :: IC_BAS =  30   ! The Baltic Sea                
  integer, parameter, public :: IC_NOS =  31   ! The North Sea                 
  integer, parameter, public :: IC_ATL =  32   ! Remaining NE Atlantic Ocean   
  integer, parameter, public :: IC_MED =  33   ! The Mediterranean Sea         
  integer, parameter, public :: IC_BLS =  34   ! The Black Sea                 
  integer, parameter, public :: IC_NAT =  35   ! Natural marine sources        
  integer, parameter, public :: IC_RUO =  36   ! Kola/Karelia                  
  integer, parameter, public :: IC_RUP =  37   ! St.Petersburg/Novgorod-Pskov  
  integer, parameter, public :: IC_RUA =  38   ! Kaliningrad                   
  integer, parameter, public ::  IC_BY =  39   ! Belarus                       
  integer, parameter, public ::  IC_UA =  40   ! Ukraine                       
  integer, parameter, public ::  IC_MD =  41   ! Moldova,                      
  integer, parameter, public :: IC_RUR =  42   ! Rest                          
  integer, parameter, public ::  IC_EE =  43   ! Estonia                       
  integer, parameter, public ::  IC_LV =  44   ! Latvia
  integer, parameter, public ::  IC_LT =  45   ! Lithuania                     
  integer, parameter, public ::  IC_CZ =  46   ! Czech                         
  integer, parameter, public ::  IC_SK =  47   ! Slovakia                      
  integer, parameter, public ::  IC_SI =  48   ! Slovenia                      
  integer, parameter, public ::  IC_HR =  49   ! Croatia                       
  integer, parameter, public ::  IC_BA =  50   ! Bosnia                        
  integer, parameter, public ::  IC_YU =  51   ! Yugoslavia                    
  integer, parameter, public ::  IC_MK =  52   ! Macedonia,                    
  integer, parameter, public ::  IC_KZ =  53   ! Kazakstan                     
  integer, parameter, public ::  IC_GE =  54   ! Georgia                       
  integer, parameter, public ::  IC_CY =  55   ! Cyprus                        
  integer, parameter, public ::  IC_AM =  56   ! Armenia                       
  integer, parameter, public ::  IC_MT =  57   ! Malta                         
  integer, parameter, public :: IC_ASI =  58   ! Other                         
  integer, parameter, public ::  IC_LI =  59   ! Lihtenstein                   
  integer, parameter, public ::  IC_DE =  60   ! Germany                       
  integer, parameter, public ::  IC_RU =  61   ! Russian                       
  integer, parameter, public ::  IC_MC =  62   ! Monaco                        
  integer, parameter, public :: IC_NOA =  63   ! North Africa                  
  integer, parameter, public ::  IC_EU =  64   ! European
  integer, parameter, public ::  IC_US =  65   ! USA
  integer, parameter, public ::  IC_CA =  66   ! Canada
  integer, parameter, public ::  IC_DUMMY1 =  67 ! Not-defined
  integer, parameter, public :: IC_KG  =  68 ! Kyrgyzstan(outside dommain)
  integer, parameter, public :: IC_AZ  =  69   ! Azerbaijan                 
  integer, parameter, public :: IC_ATX  =  70   ! ATL outside emep
  integer, parameter, public :: IC_RUX  =  71   ! RU outside emep





  contains
  !
  subroutine Country_Init()

  ! Set the country details. Note that time-zones for some areas are either
  ! difficult (Russia should be 3 to 12) or not relevant (sea areas,
  ! volcanoes). This needs to be thought about in using these figures.

  !------------------  code   index sea timezone  Name  ----------------------!
  
Country( IC_AL ) = cc(  "AL " ,  1 ,F,  1  , "Albania                       " )
Country( IC_AT ) = cc(  "AT " ,  2 ,F,  1  , "Austria                       " )
Country( IC_BE ) = cc(  "BE " ,  3 ,F,  1  , "Belgium                       " )
Country( IC_BG ) = cc(  "BG " ,  4 ,F,  2  , "Bulgaria                      " )
Country( IC_CS ) = cc(  "CS " ,  5 ,F,  1  , "Former Czechoslovakia         " )
Country( IC_DK ) = cc(  "DK " ,  6 ,F,  1  , "Denmark                       " )
Country( IC_FI ) = cc(  "FI " ,  7 ,F,  2  , "Finland                       " )
Country( IC_FR ) = cc(  "FR " ,  8 ,F,  1  , "France                        " )
Country( IC_GDR) = cc(  "GDR" ,  9 ,F,  1  , "Former East Germany           " )
Country( IC_FRG) = cc(  "FRG" , 10 ,F,  1  , "Former Fed. Rep. of Germany   " )
Country( IC_GR ) = cc(  "GR " , 11 ,F,  2  , "Greece                        " )
Country( IC_HU ) = cc(  "HU " , 12 ,F,  1  , "Hungary                       " )
Country( IC_IS ) = cc(  "IS " , 13 ,F,  0  , "Iceland                       " )
Country( IC_IE ) = cc(  "IE " , 14 ,F,  0  , "Ireland                       " )
Country( IC_IT ) = cc(  "IT " , 15 ,F,  1  , "Italy                         " )
Country( IC_LU ) = cc(  "LU " , 16 ,F,  1  , "Luxembourg                    " )
Country( IC_NL ) = cc(  "NL " , 17 ,F,  1  , "Netherlands                   " )
Country( IC_NO ) = cc(  "NO " , 18 ,F,  1  , "Norway                        " )
Country( IC_PL ) = cc(  "PL " , 19 ,F,  1  , "Poland                        " )
Country( IC_PT ) = cc(  "PT " , 20 ,F,  1  , "Portugal                      " )
Country( IC_RO ) = cc(  "RO " , 21 ,F,  2  , "Romania                       " )
Country( IC_ES ) = cc(  "ES " , 22 ,F,  1  , "Spain                         " )
Country( IC_SE ) = cc(  "SE " , 23 ,F,  1  , "Sweden                        " )
Country( IC_CH ) = cc(  "CH " , 24 ,F,  1  , "Switzerland                   " )
Country( IC_TR ) = cc(  "TR " , 25 ,F,  2  , "Turkey                        " )
Country( IC_SU ) = cc(  "SU " , 26 ,F,  3  , "Former USSR                   " )
Country( IC_GB ) = cc(  "GB " , 27 ,F,  0  , "United Kingdom                " )
Country( IC_VUL) = cc(  "VUL" , 28 ,F,  1  , "Volcanoes                     " )
Country( IC_REM) = cc(  "REM" , 29 ,F,  1  , "Remaining land areas          " )
Country( IC_BAS) = cc(  "BAS" , 30 ,T,  1  , "The Baltic Sea                " )
Country( IC_NOS) = cc(  "NOS" , 31 ,T,  1  , "The North Sea                 " )
Country( IC_ATL) = cc(  "ATL" , 32 ,T,  1  , "Remaining NE Atlantic Ocean   " )
Country( IC_MED) = cc(  "MED" , 33 ,T,  1  , "The Mediterranean Sea         " )
Country( IC_BLS) = cc(  "BLS" , 34 ,T,  1  , "The Black Sea                 " )
Country( IC_NAT) = cc(  "NAT" , 35 ,F,  1  , "Natural marine sources        " )
Country( IC_RUO) = cc(  "RUO" , 36 ,F,  3  , "Kola/Karelia                  " )
Country( IC_RUP) = cc(  "RUP" , 37 ,F,  3  , "St.Petersburg/Novgorod-Pskov  " )
Country( IC_RUA) = cc(  "RUA" , 38 ,F,  3  , "Kaliningrad                   " )
Country( IC_BY ) = cc(  "BY " , 39 ,F,  2  , "Belarus                       " )
Country( IC_UA ) = cc(  "UA " , 40 ,F,  2  , "Ukraine                       " )
Country( IC_MD ) = cc(  "MD " , 41 ,F,  2  , "Moldova, Republic of          " )
Country( IC_RUR) = cc(  "RUR" , 42 ,F,  4  , "Rest of Russia                " )
Country( IC_EE ) = cc(  "EE " , 43 ,F,  2  , "Estonia                       " )
Country( IC_LV ) = cc(  "LV " , 44 ,F,  2  , "Latvia                        " )
Country( IC_LT ) = cc(  "LT " , 45 ,F,  2  , "Lithuania                     " )
Country( IC_CZ ) = cc(  "CZ " , 46 ,F,  1  , "Czech                         " )
Country( IC_SK ) = cc(  "SK " , 47 ,F,  1  , "Slovakia                      " )
Country( IC_SI ) = cc(  "SI " , 48 ,F,  1  , "Slovenia                      " )
Country( IC_HR ) = cc(  "HR " , 49 ,F,  1  , "Croatia                       " )
Country( IC_BA ) = cc(  "BA " , 50 ,F,  1  , "Bosnia and Herzegovina        " )
Country( IC_YU ) = cc(  "YU " , 51 ,F,  1  , "Yugoslavia                    " )
Country( IC_MK ) = cc(  "MK " , 52 ,F,  1  , "Macedonia, The F.Yugo.Rep. of " )
Country( IC_KZ ) = cc(  "KZ " , 53 ,F,  6  , "Kazakstan                     " )
Country( IC_GE ) = cc(  "GE " , 54 ,F,  4  , "Georgia                       " )
Country( IC_CY ) = cc(  "CY " , 55 ,F,  2  , "Cyprus                        " )
Country( IC_AM ) = cc(  "AM " , 56 ,F,  4  , "Armenia                       " )
Country( IC_MT ) = cc(  "MT " , 57 ,F,  1  , "Malta                         " )
Country( IC_ASI) = cc(  "ASI" , 58 ,F,  0  , "Other Atlantic areas          " )
Country( IC_LI ) = cc(  "LI " , 59 ,F,  1  , "Lichtenstein                  " )
Country( IC_DE ) = cc(  "DE " , 60 ,F,  1  , "Germany                       " )
Country( IC_RU ) = cc(  "RU " , 61 ,F,  3  , "Russian Federation            " )
Country( IC_MC ) = cc(  "MC " , 62 ,F,  1  , "Monaco                        " )
Country( IC_NOA) = cc(  "NOA" , 63 ,F,  1  , "North Africa                  " )
Country( IC_EU ) = cc(  "EU " , 64 ,F,  1  , "European Community            " )
Country( IC_US ) = cc(  "US " , 65 ,F,  1  , "USA                           " )
Country( IC_CA ) = cc(  "CA " , 66 ,F,  1  , "Canada                        " )
Country( IC_DUMMY1 ) &
                 = cc(  "N/A" , 67 ,F,  -100  , "Not_defined                   " )
Country( IC_KG ) = cc(  "KG " , 68 ,F,  6  , "Kyrgyzstan                    " )
Country( IC_AZ ) = cc(  "AZ " , 69 ,F,  3  , "Azerbaijan                    " )
Country( IC_ATX )= cc(  "ATX" , 70 ,T,  1  , "Atlantic outside emep        " )
Country( IC_RUX )= cc(  "RUX" , 71 ,F,  4  , "Russian Fed. outside emep     " )



  end subroutine Country_Init
end module Country_ml
