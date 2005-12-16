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

  integer, parameter, public :: NLAND = 193
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
  integer, parameter, public ::  IC_TR =  25   ! Turkey          ME     
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
  integer, parameter, public ::  IC_KZ =  53   ! Kazakstan FE                
  integer, parameter, public ::  IC_GE =  54   ! Georgia                       
  integer, parameter, public ::  IC_CY =  55   ! Cyprus                        
  integer, parameter, public ::  IC_AM =  56   ! Armenia                       
  integer, parameter, public ::  IC_MT =  57   ! Malta                         
  integer, parameter, public :: IC_ASI =  58   ! Other                         
  integer, parameter, public ::  IC_LI =  59   ! Lichtenstein                  
  integer, parameter, public ::  IC_DE =  60   ! Germany                       
  integer, parameter, public ::  IC_RU =  61   ! Russian    FE               
  integer, parameter, public ::  IC_MC =  62   ! Monaco                        
  integer, parameter, public :: IC_NOA =  63   ! North Africa                  
  integer, parameter, public ::  IC_EU =  64   ! European
  integer, parameter, public ::  IC_US =  65   ! USA
  integer, parameter, public ::  IC_CA =  66   ! Canada
  integer, parameter, public ::  IC_DUMMY1 =  67 ! Not-defined
  integer, parameter, public :: IC_KG  =  68 ! Kyrgyzstan(outside dommain) FE
  integer, parameter, public :: IC_AZ  =  69   ! Azerbaijan                 
  ! Continent Asia
  integer, parameter, public ::  IC_AF =  70   ! Afghanistan FE
  integer, parameter, public ::  IC_BD =  71   ! Bangladesh FE
  integer, parameter, public ::  IC_BN =  72   ! Brunei FE
  integer, parameter, public ::  IC_BT =  73   ! Bhutan FE
  integer, parameter, public ::  IC_KH =  74   ! Cambodia FE
  integer, parameter, public ::  IC_LK =  75   ! Sri Lanka FE
  integer, parameter, public ::  IC_CN =  76   ! China FE
  integer, parameter, public ::  IC_MP =  77   ! Northern Mariana Islands FE
  integer, parameter, public :: IC_GZX =  78   ! Gaza Strip, The ME
  integer, parameter, public ::  IC_HK =  79   ! Hong Kong FE
  integer, parameter, public ::  IC_ID =  80   ! Indonesia FE
  integer, parameter, public ::  IC_IN =  81   ! India FE
  integer, parameter, public ::  IC_IR =  82   ! Iran ME
  integer, parameter, public ::  IC_IL =  83   ! Israel ME
  integer, parameter, public :: IC_IQX =  84   ! Iraq-saudi Arabia Neutral Zone ME
  integer, parameter, public ::  IC_IQ =  85   ! Iraq ME
  integer, parameter, public ::  IC_JP =  86   ! Japan FE
  integer, parameter, public ::  IC_JO =  87   ! Jordan ME
  integer, parameter, public ::  IC_KP =  88   ! Korea, Peoples Republic Of FE
  integer, parameter, public ::  IC_KR =  89   ! Korea, Republic Of FE
  integer, parameter, public ::  IC_KW =  90   ! Kuwait ME
  integer, parameter, public ::  IC_LA =  91   ! Laos FE
  integer, parameter, public ::  IC_LB =  92   ! Lebanon ME
  integer, parameter, public ::  IC_MN =  93   ! Mongolia FE
  integer, parameter, public ::  IC_MM =  94   ! Myanmar FE
  integer, parameter, public ::  IC_MV =  95   ! Maldives FE
  integer, parameter, public ::  IC_OM =  96   ! Oman ME
  integer, parameter, public ::  IC_MY =  97   ! Malaysia FE
  integer, parameter, public ::  IC_NP =  98   ! Nepal FE
  integer, parameter, public ::  IC_PK =  99   ! Pakistan FE
  integer, parameter, public ::  IC_QA = 100   ! Qatar ME
  integer, parameter, public ::  IC_PH = 101   ! Philippines, The FE
  integer, parameter, public ::  IC_SA = 102   ! Saudi Arabia ME
  integer, parameter, public ::  IC_SG = 103   ! Singapore FE
  integer, parameter, public ::  IC_SY = 104   ! Syria ME
  integer, parameter, public ::  IC_TJ = 105   ! Tajikistan FE
  integer, parameter, public ::  IC_AE = 106   ! United Arab Emirates, The ME
  integer, parameter, public ::  IC_TH = 107   ! Thailand FE
  integer, parameter, public ::  IC_TM = 108   ! Turkmenistan FE
  integer, parameter, public ::  IC_TW = 109   ! Taiwan FE
  integer, parameter, public ::  IC_UZ = 110   ! Uzbekistan FE
  integer, parameter, public ::  IC_VN = 111   ! Vietnam FE
  integer, parameter, public ::  IC_YE = 112   ! Yemen ME
  ! Continent Africa
  integer, parameter, public ::  IC_DZ = 113   ! Algeria
  integer, parameter, public ::  IC_AO = 114   ! Angola
  integer, parameter, public ::  IC_BJ = 115   ! Benin
  integer, parameter, public ::  IC_BI = 116   ! Burundi
  integer, parameter, public ::  IC_TD = 117   ! Chad
  integer, parameter, public ::  IC_CG = 118   ! Congo, The
  integer, parameter, public ::  IC_CD = 119   ! Zaire
  integer, parameter, public ::  IC_CM = 120   ! Cameroon
  integer, parameter, public ::  IC_CF = 121   ! Central African Republic
  integer, parameter, public ::  IC_CV = 122   ! Cape Verde
  integer, parameter, public ::  IC_DJ = 123   ! Djibouti
  integer, parameter, public ::  IC_EG = 124   ! Egypt
  integer, parameter, public ::  IC_GQ = 125   ! Equatorial Guinea
  integer, parameter, public ::  IC_ET = 126   ! Ethiopia
  integer, parameter, public ::  IC_GM = 127   ! Gambia, The
  integer, parameter, public ::  IC_GA = 128   ! Gabon
  integer, parameter, public ::  IC_GH = 129   ! Ghana
  integer, parameter, public ::  IC_GN = 130   ! Guinea
  integer, parameter, public ::  IC_CI = 131   ! Ivory Coast, The
  integer, parameter, public ::  IC_KE = 132   ! Kenya
  integer, parameter, public ::  IC_LR = 133   ! Liberia
  integer, parameter, public ::  IC_LY = 134   ! Libya
  integer, parameter, public ::  IC_ML = 135   ! Mali
  integer, parameter, public ::  IC_MA = 136   ! Morocco
  integer, parameter, public ::  IC_MR = 137   ! Mauritania
  integer, parameter, public ::  IC_NE = 138   ! Niger
  integer, parameter, public ::  IC_NG = 139   ! Nigeria
  integer, parameter, public ::  IC_GW = 140   ! Guinea-bissau
  integer, parameter, public ::  IC_RW = 141   ! Rwanda
  integer, parameter, public ::  IC_SN = 142   ! Senegal
  integer, parameter, public ::  IC_SL = 143   ! Sierra Leone
  integer, parameter, public ::  IC_ST = 144   ! Sao Tome And Principe
  integer, parameter, public ::  IC_SO = 145   ! Somalia
  integer, parameter, public ::  IC_SD = 146   ! Sudan
  integer, parameter, public ::  IC_TG = 147   ! Togo
  integer, parameter, public ::  IC_TN = 148   ! Tunisia
  integer, parameter, public ::  IC_TZ = 149   ! Tanzania
  integer, parameter, public ::  IC_UG = 150   ! Uganda
  integer, parameter, public ::  IC_BF = 151   ! Burkina Faso
  integer, parameter, public ::  IC_NA = 152   ! Namibia
  integer, parameter, public ::  IC_EH = 153   ! Western Sahara
  integer, parameter, public ::  IC_ZM = 154   ! Zambia
  ! Continent North America
  integer, parameter, public ::  IC_AG = 155   ! Antigua And Barbuda
  integer, parameter, public ::  IC_BB = 156   ! Barbados
  integer, parameter, public ::  IC_BS = 157   ! Bahamas, The
  integer, parameter, public ::  IC_BZ = 158   ! Belize
  integer, parameter, public ::  IC_BM = 159   ! Bermuda
  integer, parameter, public ::  IC_CR = 160   ! Costa Rica
  integer, parameter, public ::  IC_CU = 161   ! Cuba
  integer, parameter, public ::  IC_DM = 162   ! Dominica
  integer, parameter, public ::  IC_DO = 163   ! Dominican Republic
  integer, parameter, public ::  IC_SV = 164   ! El Salvador
  integer, parameter, public ::  IC_GL = 165   ! Greenland
  integer, parameter, public ::  IC_GD = 166   ! Grenada
  integer, parameter, public ::  IC_GP = 167   ! Guadeloupe
  integer, parameter, public ::  IC_GT = 168   ! Guatemala
  integer, parameter, public ::  IC_HT = 169   ! Haiti
  integer, parameter, public ::  IC_HN = 170   ! Honduras
  integer, parameter, public ::  IC_JM = 171   ! Jamaica
  integer, parameter, public ::  IC_MQ = 172   ! Martinique
  integer, parameter, public ::  IC_MX = 173   ! Mexico
  integer, parameter, public ::  IC_NI = 174   ! Nicaragua
  integer, parameter, public ::  IC_PA = 175   ! Panama
  integer, parameter, public ::  IC_PR = 176   ! Puerto Rico
  integer, parameter, public ::  IC_KN = 177   ! St. Christopher-nevis
  integer, parameter, public ::  IC_VC = 178   ! St. Vincent And The Grenadines
  integer, parameter, public ::  IC_LC = 179   ! St. Lucia
  integer, parameter, public ::  IC_TT = 180   ! Trinidad And Tobago
  integer, parameter, public ::  IC_TC = 181   ! Turks And Caicos Islands
  ! Continent South America
  integer, parameter, public ::  IC_BO = 182   ! Bolivia
  integer, parameter, public ::  IC_BR = 183   ! Brazil
  integer, parameter, public ::  IC_CO = 184   ! Colombia
  integer, parameter, public ::  IC_EC = 185   ! Ecuador
  integer, parameter, public ::  IC_GF = 186   ! French Guiana
  integer, parameter, public ::  IC_GY = 187   ! Guyana
  integer, parameter, public ::  IC_SR = 188   ! Surinam
  integer, parameter, public ::  IC_PE = 189   ! Peru
  integer, parameter, public ::  IC_VE = 190   ! Venezuela
  ! Continent Australia
  integer, parameter, public ::  IC_KI = 191   ! Kiribati
  integer, parameter, public ::  IC_WS = 192   ! Western Samoa
  ! Continent non
  integer, parameter, public :: IC_NXX = 193   ! Non




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
!Country( IC_RUR) = cc(  "RUR" , 42 ,F,  4  , "Rest of Russia                " )
Country( IC_RUR) = cc(  "RUR" , 42 ,F,-100 , "Rest of Russia                " )
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
                 = cc(  "N/A" , 67 ,F,  0  , "Not_defined                   " )
Country( IC_KG ) = cc(  "KG " , 68 ,F,  6  , "Kyrgyzstan                    " )
Country( IC_AZ ) = cc(  "AZ " , 69 ,F,  3  , "Azerbaijan                    " )
  ! Continent Asia
Country( IC_AF ) = cc(  "AF ",11001,F,-100 , "Afghanistan                   " )
Country( IC_BD ) = cc(  "BD ",11002,F,-100 , "Bangladesh                    " )
Country( IC_BN ) = cc(  "BN ",11003,F,-100 , "Brunei                        " )
Country( IC_BT ) = cc(  "BT ",11005,F,-100 , "Bhutan                        " )
Country( IC_KH ) = cc(  "KH ",11006,F,-100 , "Cambodia                      " )
Country( IC_LK ) = cc(  "LK ",11007,F,-100 , "Sri Lanka                     " )
Country( IC_CN ) = cc(  "CN ",11008,F,-100 , "China                         " )
Country( IC_MP ) = cc(  "MP ",11009,F,-100 , "Northern Mariana Islands      " )
Country( IC_GZX) = cc(  "GZX",11010,F,-100 , "Gaza Strip, The               " )
Country( IC_HK ) = cc(  "HK ",11011,F,-100 , "Hong Kong                     " )
Country( IC_ID ) = cc(  "ID ",11012,F,-100 , "Indonesia                     " )
Country( IC_IN ) = cc(  "IN ",11013,F,-100 , "India                         " )
Country( IC_IR ) = cc(  "IR ",11014,F,-100 , "Iran                          " )
Country( IC_IL ) = cc(  "IL ",11015,F,-100 , "Israel                        " )
Country( IC_IQX) = cc(  "IQX",11016,F,-100 , "Iraq-saudi Arabia Neutral Zone" )
Country( IC_IQ ) = cc(  "IQ ",11017,F,-100 , "Iraq                          " )
Country( IC_JP ) = cc(  "JP ",11018,F,-100 , "Japan                         " )
Country( IC_JO ) = cc(  "JO ",11019,F,-100 , "Jordan                        " )
Country( IC_KP ) = cc(  "KP ",11020,F,-100 , "Korea, Peoples Republic Of    " )
Country( IC_KR ) = cc(  "KR ",11021,F,-100 , "Korea, Republic Of            " )
Country( IC_KW ) = cc(  "KW ",11022,F,-100 , "Kuwait                        " )
Country( IC_LA ) = cc(  "LA ",11023,F,-100 , "Laos                          " )
Country( IC_LB ) = cc(  "LB ",11024,F,-100 , "Lebanon                       " )
Country( IC_MN ) = cc(  "MN ",11025,F,-100 , "Mongolia                      " )
Country( IC_MM ) = cc(  "MM ",11026,F,-100 , "Myanmar                       " )
Country( IC_MV ) = cc(  "MV ",11027,F,-100 , "Maldives                      " )
Country( IC_OM ) = cc(  "OM ",11028,F,-100 , "Oman                          " )
Country( IC_MY ) = cc(  "MY ",11029,F,-100 , "Malaysia                      " )
Country( IC_NP ) = cc(  "NP ",11030,F,-100 , "Nepal                         " )
Country( IC_PK ) = cc(  "PK ",11031,F,-100 , "Pakistan                      " )
Country( IC_QA ) = cc(  "QA ",11032,F,-100 , "Qatar                         " )
Country( IC_PH ) = cc(  "PH ",11033,F,-100 , "Philippines, The              " )
Country( IC_SA ) = cc(  "SA ",11034,F,-100 , "Saudi Arabia                  " )
Country( IC_SG ) = cc(  "SG ",11035,F,-100 , "Singapore                     " )
Country( IC_SY ) = cc(  "SY ",11036,F,-100 , "Syria                         " )
Country( IC_TJ ) = cc(  "TJ ",11037,F,-100 , "Tajikistan                    " )
Country( IC_AE ) = cc(  "AE ",11038,F,-100 , "United Arab Emirates, The     " )
Country( IC_TH ) = cc(  "TH ",11039,F,-100 , "Thailand                      " )
Country( IC_TM ) = cc(  "TM ",11040,F,-100 , "Turkmenistan                  " )
Country( IC_TW ) = cc(  "TW ",11041,F,-100 , "Taiwan                        " )
Country( IC_UZ ) = cc(  "UZ ",11042,F,-100 , "Uzbekistan                    " )
Country( IC_VN ) = cc(  "VN ",11043,F,-100 , "Vietnam                       " )
Country( IC_YE ) = cc(  "YE ",11044,F,-100 , "Yemen                         " )
  ! Continent Africa
Country( IC_DZ ) = cc(  "DZ ",12001,F,-100 , "Algeria                       " )
Country( IC_AO ) = cc(  "AO ",12002,F,-100 , "Angola                        " )
Country( IC_BJ ) = cc(  "BJ ",12003,F,-100 , "Benin                         " )
Country( IC_BI ) = cc(  "BI ",12004,F,-100 , "Burundi                       " )
Country( IC_TD ) = cc(  "TD ",12005,F,-100 , "Chad                          " )
Country( IC_CG ) = cc(  "CG ",12006,F,-100 , "Congo, The                    " )
Country( IC_CD ) = cc(  "CD ",12007,F,-100 , "Zaire                         " )
Country( IC_CM ) = cc(  "CM ",12008,F,-100 , "Cameroon                      " )
Country( IC_CF ) = cc(  "CF ",12009,F,-100 , "Central African Republic      " )
Country( IC_CV ) = cc(  "CV ",12010,F,-100 , "Cape Verde                    " )
Country( IC_DJ ) = cc(  "DJ ",12011,F,-100 , "Djibouti                      " )
Country( IC_EG ) = cc(  "EG ",12012,F,-100 , "Egypt                         " )
Country( IC_GQ ) = cc(  "GQ ",12013,F,-100 , "Equatorial Guinea             " )
Country( IC_ET ) = cc(  "ET ",12014,F,-100 , "Ethiopia                      " )
Country( IC_GM ) = cc(  "GM ",12015,F,-100 , "Gambia, The                   " )
Country( IC_GA ) = cc(  "GA ",12016,F,-100 , "Gabon                         " )
Country( IC_GH ) = cc(  "GH ",12017,F,-100 , "Ghana                         " )
Country( IC_GN ) = cc(  "GN ",12018,F,-100 , "Guinea                        " )
Country( IC_CI ) = cc(  "CI ",12019,F,-100 , "Ivory Coast, The              " )
Country( IC_KE ) = cc(  "KE ",12020,F,-100 , "Kenya                         " )
Country( IC_LR ) = cc(  "LR ",12021,F,-100 , "Liberia                       " )
Country( IC_LY ) = cc(  "LY ",12022,F,-100 , "Libya                         " )
Country( IC_ML ) = cc(  "ML ",12023,F,-100 , "Mali                          " )
Country( IC_MA ) = cc(  "MA ",12024,F,-100 , "Morocco                       " )
Country( IC_MR ) = cc(  "MR ",12025,F,-100 , "Mauritania                    " )
Country( IC_NE ) = cc(  "NE ",12026,F,-100 , "Niger                         " )
Country( IC_NG ) = cc(  "NG ",12027,F,-100 , "Nigeria                       " )
Country( IC_GW ) = cc(  "GW ",12028,F,-100 , "Guinea-bissau                 " )
Country( IC_RW ) = cc(  "RW ",12029,F,-100 , "Rwanda                        " )
Country( IC_SN ) = cc(  "SN ",12030,F,-100 , "Senegal                       " )
Country( IC_SL ) = cc(  "SL ",12031,F,-100 , "Sierra Leone                  " )
Country( IC_ST ) = cc(  "ST ",12032,F,-100 , "Sao Tome And Principe         " )
Country( IC_SO ) = cc(  "SO ",12033,F,-100 , "Somalia                       " )
Country( IC_SD ) = cc(  "SD ",12034,F,-100 , "Sudan                         " )
Country( IC_TG ) = cc(  "TG ",12035,F,-100 , "Togo                          " )
Country( IC_TN ) = cc(  "TN ",12036,F,-100 , "Tunisia                       " )
Country( IC_TZ ) = cc(  "TZ ",12037,F,-100 , "Tanzania                      " )
Country( IC_UG ) = cc(  "UG ",12038,F,-100 , "Uganda                        " )
Country( IC_BF ) = cc(  "BF ",12039,F,-100 , "Burkina Faso                  " )
Country( IC_NA ) = cc(  "NA ",12040,F,-100 , "Namibia                       " )
Country( IC_EH ) = cc(  "EH ",12041,F,-100 , "Western Sahara                " )
Country( IC_ZM ) = cc(  "ZM ",12042,F,-100 , "Zambia                        " )
  ! Continent North America
Country( IC_AG ) = cc(  "AG ",13001,F,-100 , "Antigua And Barbuda           " )
Country( IC_BB ) = cc(  "BB ",13002,F,-100 , "Barbados                      " )
Country( IC_BS ) = cc(  "BS ",13003,F,-100 , "Bahamas, The                  " )
Country( IC_BZ ) = cc(  "BZ ",13004,F,-100 , "Belize                        " )
Country( IC_BM ) = cc(  "BM ",13005,F,-100 , "Bermuda                       " )
Country( IC_CR ) = cc(  "CR ",13006,F,-100 , "Costa Rica                    " )
Country( IC_CU ) = cc(  "CU ",13007,F,-100 , "Cuba                          " )
Country( IC_DM ) = cc(  "DM ",13008,F,-100 , "Dominica                      " )
Country( IC_DO ) = cc(  "DO ",13009,F,-100 , "Dominican Republic            " )
Country( IC_SV ) = cc(  "SV ",13010,F,-100 , "El Salvador                   " )
Country( IC_GL ) = cc(  "GL ",  601,F,-100 , "Greenland                     " )
Country( IC_GD ) = cc(  "GD ",13011,F,-100 , "Grenada                       " )
Country( IC_GP ) = cc(  "GP ",13012,F,-100 , "Guadeloupe                    " )
Country( IC_GT ) = cc(  "GT ",13013,F,-100 , "Guatemala                     " )
Country( IC_HT ) = cc(  "HT ",13014,F,-100 , "Haiti                         " )
Country( IC_HN ) = cc(  "HN ",13015,F,-100 , "Honduras                      " )
Country( IC_JM ) = cc(  "JM ",13016,F,-100 , "Jamaica                       " )
Country( IC_MQ ) = cc(  "MQ ",13017,F,-100 , "Martinique                    " )
Country( IC_MX ) = cc(  "MX ",13018,F,-100 , "Mexico                        " )
Country( IC_NI ) = cc(  "NI ",13019,F,-100 , "Nicaragua                     " )
Country( IC_PA ) = cc(  "PA ",13020,F,-100 , "Panama                        " )
Country( IC_PR ) = cc(  "PR ",13021,F,-100 , "Puerto Rico                   " )
Country( IC_KN ) = cc(  "KN ",13022,F,-100 , "St. Christopher-nevis         " )
Country( IC_VC ) = cc(  "VC ",13023,F,-100 , "St. Vincent And The Grenadines" )
Country( IC_LC ) = cc(  "LC ",13024,F,-100 , "St. Lucia                     " )
Country( IC_TT ) = cc(  "TT ",13025,F,-100 , "Trinidad And Tobago           " )
Country( IC_TC ) = cc(  "TC ",13026,F,-100 , "Turks And Caicos Islands      " )
  ! Continent South America
Country( IC_BO ) = cc(  "BO ",14001,F,-100 , "Bolivia                       " )
Country( IC_BR ) = cc(  "BR ",14002,F,-100 , "Brazil                        " )
Country( IC_CO ) = cc(  "CO ",14003,F,-100 , "Colombia                      " )
Country( IC_EC ) = cc(  "EC ",14004,F,-100 , "Ecuador                       " )
Country( IC_GF ) = cc(  "GF ",14005,F,-100 , "French Guiana                 " )
Country( IC_GY ) = cc(  "GY ",14006,F,-100 , "Guyana                        " )
Country( IC_SR ) = cc(  "SR ",14007,F,-100 , "Surinam                       " )
Country( IC_PE ) = cc(  "PE ",14008,F,-100 , "Peru                          " )
Country( IC_VE ) = cc(  "VE ",14009,F,-100 , "Venezuela                     " )
  ! Continent Australia
Country( IC_KI ) = cc(  "KI ",15001,F,-100 , "Kiribati                      " )
Country( IC_WS ) = cc(  "WS ",15002,F,-100 , "Western Samoa                 " )
  ! Continent non
Country( IC_NXX) = cc(  "NXX", 9999,F,-100 , "Non                           " )

  end subroutine Country_Init
end module Country_ml
