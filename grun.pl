#!/usr/local/bin/perl
######################################################################
# Features
# 1. work directory now deduced from user-name
# 2. Check that the number of processors asked for by bsub is same as
#   number set here. Note that this method will not allow a simple
#   interactive use of grun.pl, except if COMPILE_ONLY is set to one.
#   "bsub" must be used.
# 3. Produces Remove.sh file which contains "rm xxxx" for all files
#    which should have been deleted at the end of the run. If a run
#    hangs or crashes for some reason, type "sh Remove.sh" to get rid
#    of all the linked files.  
# Bugs
# None known?
######################################################################
#
# Script to prepare files before running Eulerian model on 
# Odin. The main advantage of this script is that
# the domain size  and input/output files and directories
# can be easily changed. The script does this by modifying the
# Makefile.pat (on huge only) and Par_ml.pat file to produce the corresponding
# Makefile and Par_ml.f90 file with appropriate set-up.  
#
# To submit on GRIDUR, put grun.pl into work directory, and from there
# do something like:
#
#    bsub -n nr_of_nodes -o log.out < grun.pl

#
#   - assigns felt files to filxxx in the current directory
#
# It is assumed that all the sequential felt files have names on the
# form f$hh.$yy$mm$dd
#
######################################################################

#------------------------------------------------------------------------------
# perl version of eulmod driver script
# (C) 1999 Jorn Amundsen, NTNU
# modified to run unified models, D. Simpson, 2000
# mylink added and clean-up, ds, 2001
#------------------------------------------------------------------------------

system "/usr/bin/ja"; #start of resource info gathering
require "flush.pl";

# <---------- "Pattern files - different on huge and gridur ---->
# These are changed by the  script accroding to the domain and
# number of processors selected.
# (but note that this script file is onlt for gridur these days)
#huge: @pattern_files = ( Makefile.pat Par_ml.pat );

@pattern_files = qw ( Par_ml.pat );

@bigfile_list = qw ( ); #list of files to be compressed

# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables 
#      are explained below, and derived variables set later.-

$year = "2003";
( $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

my $SR = 0;   #NEW Set to 1 for source-receptor stuff
my $PM_ADDED     = 1;  # Adds in PM emissions from NOx inventory scaling
my $AFRICA_ADDED = 1;  # Adds in African emissions for y=1..11
my $MERLIN_CITY= 0;  # Adds in African emissions for y=1..11

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

$iyr_trend = $year;  
$iyr_trend = "2010" if $SR ;  # 2010 assumed for SR runs here

print "Year is $yy YEAR $year Trend year $ir_trend\n";

if ( $year == 2000 ) {
  $MetDir = "/work/mifads/metdata/$year" ;
} elsif ( $year == 1997 ) {
  $MetDir = "/work/mifapw/metdata/$year" ;
} elsif ( $year == 2002 ) {
  $MetDir = "/work/mifahf/metdata/$year" ;
} elsif ( $year == 1995 ) {
  $MetDir = "/work/mifads/metdata/$year" ;  # Needed for 1997
} else { # 2003 etc.
  $MetDir = "/work/mifapw/metdata/$year" ;  # Default
}

#---  User-specific directories (changeable)

$DAVE        = "/home/u2/mifads";      
$JOFFEN      = "/home/u2/mifajej";      
$HILDE       = "/home/u6/mifahf";      
$STEFFEN     = "/home/u2/mifaung";      
$SVETLANA    = "/home/u2/mifast";      
$PETER       = "/home/u4/mifapw";      

$USER        =  $PETER ;      

$USER  =~ /(\w+ $)/x ;       # puts word characters (\w+) at end ($) into "$1"
$WORK{$USER} = "/work/$1";   # gives e.g. /work/mifads


#ds - simplified treatment of BCs and emissions:

$OZONE = 1, $ACID = 0;     # Specify model type here

  # check:
  die "Must choose ACID or OZONE" if ( $OZONE+$ACID>1 or $OZONE+$ACID==0 );



# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

if ( $OZONE ) {
    #$OZONEDIR    = "$HILDE/BC_data/CTM2_data/50Data"; 
     $OZONEDIR    = "$HILDE/BC_data/LOGAN_O3_DATA/50Data_900mbar"; 
    #$OZONEDIR    = "$HILDE/BC_data/Fortuin_data/50Data"; 
     @emislist = qw ( sox nox nh3 co voc pm25 pmco ); 
     $testv       = "rv2_0_5";
     $runlabel1    = "$testv";           # NO SPACES! SHORT name (used in CDF names)
     $runlabel2    = "${testv}_$year";   # NO SPACES! LONG (written into CDF files)
  $BCDIR    = "$HILDE/BC_data/Unimod.rv2_0beta.1997";  # for CH3COO2 and H2O2
  $BCDIRO3    = "$HILDE/BC_data/LOGAN_O3_DATA/50Data_900mbar_mixratio";  # Logan for O3

} elsif ( $ACID ) {
      $OZONEDIR    = "$HILDE/BC_data/Unimod.rv2_0beta.1997";
     @emislist = qw ( sox nox nh3 pm25 pmco ) ;
     $testv       = "rv2_0beta.oh.acid";
     $runlabel1    = "TEST_of_$testv";   # NO SPACES! SHORT name (used in CDF names)
     $runlabel2    = "${testv}_XXX_$year";   # NO SPACES! LONG (written into CDF files)
  #XD: New source of BCs for OH, CH3COO2, O3 and H2O2: (from OZONE)
  $BCDIR    = "$HILDE/BC_data/Unimod.rv2_0beta.${year}";
  $BCDIRO3    = $BCDIR;  # 
  }
  die "NO DIR $BCDIR \n" unless -d $BCDIR;
  die "NO DIR $BCDIRO3 \n" unless -d $BCDIRO3;
$version     = "Unimod" ;  
$subv        = "$testv" ;                  # sub-version (to track changes)
$Split       = "BASE_MAR2004" ;               #  -- Scenario label for MACH - DS
$ProgDir     = "$USER/Unify/Unimod.$testv";   # input of source-code
#$ProgDir     = "$USER/Unify/temp";   # input of source-code
$MyDataDir   = "$USER/Unify/MyData";          # for each user's femis, etc.
$DaveDataDir   = "$DAVE/Unify/MyData";          # for each user's femis, etc.
$DataDir     = "$DAVE/Unify/Data";      # common files, e.g. ukdep_biomass.dat
$PROGRAM     = "$ProgDir/$version";         # programme
$WORKDIR     = "$WORK{$USER}/Unimod.$testv";    # working directory

# $femis       = "$MyDataDir/femis.dat";      # emission control file
$femis_dir   = "$WORKDIR/D_femis";      # emission control file
system("mkdir $femis_dir");
$Africa      = "$DAVE/Unify/D_emis";        # Emissions for Africa, y=1..11
$timeseries  = "$DAVE/Unify/D_timeseries";   # New timeseries (ds 14/1/2003) 

#ds new for SR version. One lement arrays if SR not specified:
$mk_SRfemis = "$DAVE/Unify/Progs/mkp.SRfemis";   #  Generates femis files

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later. 
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR
@runs        = ( $scenario );
@scenarios   = ( "Base");
system("cp  $MyDataDir/femis.dat  $femis_dir/femis.$scenario");  


$emisdir     = "$PETER/Vigdis/Emissions/Modruns/Modrun03";
if ( $SR ) {
	#? $emisdir     = "$emisdir/2003_emis2010_CLE_2000_V5";    # emissions
	$emisdir     = "$emisdir/2004_emis2010_CLE_2000_V6";    # emissions
} elsif ( $year == 2000)  {
	#? $emisdir     = "$emisdir/2003_emis00_V5_WEB";    # emissions
	#$emisdir     = "$emisdir/2004_emis00_V6";    # emissions
       $emisdir     = "$SVETLANA/Unify/MyData/emission/2004_emis00_V2";    # emissions
} elsif ( $year == 2001)  {
	$emisdir     = "$emisdir/2004_emis01_V6";    # emissions
} elsif ( $year == 2002)  {
       $emisdir     = "$SVETLANA/Unify/MyData/emission/2004_emis02_V2";    # emissions
} elsif ( $year == 2003)  {
       $emisdir     = "$SVETLANA/Unify/MyData/emission/2004_emis02_V2";    # emissions
} else {
	#
       $emisdir     = "$SVETLANA/Unify/MyData/emission/2004_emis00_V2"; # 2000 emissions def
}

# Specify small domain if required. 
#                 x0   x1  y0   y1
@largedomain = (   1, 170,  1, 133 ) ;
#@smalldomain = ( 101, 140, 51,  90 ) ;      # (changeable)
#@smalldomain = (  71, 150, 31, 100 ) ;      # (changeable)
#@smalldomain = (  95, 115, 46, 66 ) ;      # ERROR search (changeable)
#@smalldomain = (  36, 160, 11, 123 ) ;      # (changeable)
#@smalldomain = (  20, 167,  1, 122 ) ;    # OSPAR/HELCOM domain
@smalldomain = (  18, 169,  7, 124 ) ;     # OSPAR/HELCOM domain+border-south
#@smalldomain = (  60, 169,  40, 124 ) ;     # OSPAR/HELCOM domain+border-south
#@smalldomain = (  36, 167, 12, 122 ) ;    # EMEP domain
#@smalldomain = (  36, 130, 31, 123 ) ;      # (changeable)
#@smalldomain = (  39, 120, 31, 123 ) ;      # (changeable)
#@smalldomain = @largedomain ;     # If you want to run for the whole domain, 
                                    # simply uncomment this 

$RESET        = 0   ;  # usually 0 (false) is ok, but set to 1 for full restart
$COMPILE_ONLY = 0   ;  # usually 0 (false) is ok, but set to 1 for compile-only
$INTERACTIVE  = 0   ;  # usually 0 (false), but set to 1 to make program stop
                       # just before execution - so code can be run interactivel.

$NDX   = 8;           # Processors in x-direction
$NDY   =  4;           # Processors in y-direction
if ( $INTERACTIVE ) { $NDX = $NDY = 1 };


@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

$mm1   =  1;       # first month
$mm2   =  1 ;       # last month
$NTERM_CALC =  calc_nterm($mm1,$mm2);

$NTERM =   $NTERM_CALC;    # sets NTERM for whole time-period
  # -- or --
 $NTERM = 9;       # for testing, simply reset here

  print "NTERM_CALC = $NTERM_CALC, Used NTERM = $NTERM\n";

# <---------- end of normal use section ---------------------->
# <---------- start of SR   use section ---------------------->

#SR additions:
# Define all coutnries and nums here: Heiko's lib might be faster:
%country_nums = ( AT => 2, BE => 3, DK => 6, FI => 7, FR => 8,
                  DE => 10, GR => 11, IE => 14, IT => 15,
                  NL => 17 , PT => 20, ES => 22, SE => 23,
                  UK => 27,
                  HU => 12, PL => 19, CY => 55, CZ => 46, EE => 43,
                  LT => 45, LV => 44, MA => 57, SK => 47, SI => 48 );
# EU countries:
@eu = qw ( AT BE DK FI FR DE GR IE IT NL PT ES SE UK );
@eu1 = qw ( AT BE DK FI FR );  # DE GR );
@eu2 = qw ( IE IT NL PT );
@eu3 = qw ( ES SE UK DE GR );  # Added 2 from eu1
@euaccess = qw ( HU PL CY CZ EE LT LV MA SK SI );
@eu25 = ( @eu, @euaccess );
@countries = @eu3 ;  # List of wanted countries this run
@countries = qw ( IT DE ) ;  # List of wanted countries this run
                                # EU25 is "special"
my $scenario = "Base";
if ( $SR ) {
        $base        = "CLE";
        $Split       = "CLE_MAR2004";    # IER VOC splits
        $rednflag    = "P15";  # 10% reduction for label
        $redn        = "0.85";  # 10% reduction
        @polls       = qw ( BASE NP );  # NP, BASE already done

	$nrun = 0;
	foreach $pollut ( @polls ) {
	    foreach $cc ( @countries ) {

       		$scenario = "${base}_${cc}_${pollut}_${rednflag}";
       		$scenario = "${base}" if  $pollut eq "BASE";

		# create list of runs to be done:
		$runs[$nrun] = $scenario;
       		print "scenario: N $nrun S $scenario \n";
		$nrun++ ;
	
		last if $pollut eq "BASE";  # Do not repeat for each cc
	    }
	}
	# Run mkp.SRfemis to generate femis files. This script will die 

	system("$mk_SRfemis -b $base -f $rednflag -r $redn -d $femis_dir");
}


# <---------- end of user-changeable section ----------------->
#               (normally, that is...)

#--- Other eulmod configs

$NSCAL =  0;
$NASS  =  0;        # Set to one if "dump" of all concentrations wanted at end
# 
# Check that we have an existing prog dir:
die "Wrong ProgDir\n" unless -d $ProgDir;

# quick check that the small domain
# is within the current 170, 133 large domain

die " -- Domain error!!!" if ( 
    $smalldomain[0] < $largedomain[0] ||  $smalldomain[1] > $largedomain[1] || 
    $smalldomain[2] < $largedomain[2] ||  $smalldomain[3] > $largedomain[3] ) ;


#--- calculate number of processors

$NPROC =  $NDX * $NDY ;

# and check that we have the same from the bsub command, by
# accessing the LSB_MCPU_HOSTS environment variable.
# (For a 5 cpu job this looks like: LSB_MCPU_HOSTS=gridur.ntnu.no 5)

$nproc_bsub = (split/\s+/,$ENV{'LSB_MCPU_HOSTS'})[1];

if ( ! $INTERACTIVE && ! $COMPILE_ONLY  && $NPROC != $nproc_bsub ) {
    die " -- Requested wrong number of processors --
      bsub asked for $nproc_bsub whereas NPROC = $NDX x $NDY = $NPROC \n";
}


# ---- calculate domain widths
# (For the model, we need first x0, y0, then width (number of cells) in x and y)

$dom_x0 = $smalldomain[0]; 
$dom_y0 = $smalldomain[2]; 
$dom_wx = $smalldomain[1] - $smalldomain[0] + 1 ; 
$dom_wy = $smalldomain[3] - $smalldomain[2] + 1 ; 


#---  Data directories   (do not change so often.... )

$MIROOT      = "/home/";
$MIFAJEJ     = "/home/u2/mifajej";     # Needed for some input data
$SRCINPUTDIR = "/home/u2/mifajej/data/climatology";   # landuse, snow, rough....
#FIX#$BCINPUTDIR  = "$HILDE/BC_data/LOGAN_O3_DATA/150Data"; #  Logan boundary conditions
$o3dir       = "$MIFAJEJ/data" ;       # ancat (=aircraft), 
$MIFADS      = "/home/u2/mifajej/data/mifads";    # Ozone data for now
$EULDATA     = "/home/u2/mifajej/data/euldata" ;   # emissions, monthlyfac,..., JAN.dat*

@MIDIRS  = ("$MIROOT", "$USER", "$MIFAJEJ", "$MIFADS" , "$o3dir", "$EULDATA");

#--- Verify data directories
foreach $d ( @MIDIRS, $SRCINPUTDIR, $WORKDIR ) {
    unless ( -d "$d" &&  -x _ && -r _ ) {
	die "*** ERROR *** directory $d not accessible. Exiting.\n";
    }
}

#--- Verify New pattern files

foreach $f ( @pattern_files ) {  
    unless ( -r "$ProgDir/$f" ) {
	die "*** ERROR *** file $f not available. Exiting.\n";
    }
}


#--- Calendar stuff
@month_abbrev = ('','jan','feb','mar','apr','may','jun','jul','aug','sep','oct',
		 'nov','dec');
@season_abbrev = ('','jan','jan','apr','apr','apr','jul','jul','jul','oct','oct',
		 'oct','jan');
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);

%seasons = ( "jan" => "01", "apr" => "02", "jul" => "03" , "oct" => "04") ;


#--- adjust for leap year
$month_days[2] += leap_year($year);

# --- Start new compilations if needed

# --- We check if the already-compiled version of $PROGRAM is the same
#     as the one we are asking for. The old configuration should have 
#     been stored in Make.log.  Read in this and split terms into 
#     : model, emis, ndx, ndy and domain sizes


chdir "$ProgDir";

#-- generate Makefile each time, to avoid forgetting changed "pat" file!
#odin system "sed 's/NPROC/$NPROC/'  Makefile.pat > Makefile" ;

if ( $RESET ) { unlink ("Make.log") }  # Gone!

open(MAKELOG,"<Make.log");
($oldversion, $olddx , $olddy , $old_x0, $old_y0, $old_wx, $old_wy ) = split ' ', <MAKELOG> ;
close(MAKELOG);

print "From Make.log we had Subversion $subv
           Procs:  $olddx $olddy 
           Domain: $old_x0, $old_y0, $old_wx, $old_wy \n";


if ( $oldversion ne $subv ) {

   print " We are changing version!!!!!............. 
            from $oldversion to $subv \n " ;
   $RESET = 1 ;
}

if ( $NDX      != $olddx  || $NDY      != $olddy  ||
     $dom_x0   != $old_x0 || $dom_y0   != $old_y0 ||
     $dom_wx   != $old_wx || $dom_wy   != $old_wy     ) {

   print "Need to re-compile for new processor or domain  setup: 
           Procs:  $NDX $NDY
           Domain: $dom_x0, $dom_y0, $dom_wx, $dom_wy \n";
   $RESET = 1 ;
}


if ( $RESET == 1  ) { ########## Recompile everything!


   # Set values for domain size in Par_ml.f90 : 
   foreach $f ( qw ( Par_ml )) {
      open(EULPAR,"<$f.pat");
      open(EULOUT,">$f.f90");
      print "changing domain, nproc in $f.pat to $f.f90 \n";
      while ($line=<EULPAR>) {
          $line =~ s/nprocx/$NDX/ ; 
          $line =~ s/nprocy/$NDY/ ;
          $line =~ s/domainx0/$dom_x0/ ;
          $line =~ s/domainy0/$dom_y0/ ;
          $line =~ s/domaindx/$dom_wx/ ;
          $line =~ s/domaindy/$dom_wy/ ;
          print EULOUT $line ;
      }
      close(EULOUT) ;
   }
   # For now, we simply recompile everything!
   unlink($PROGRAM);
   system "touch -c  *.f *.f90 *.F *.F90";
}
system "pwd";
print "Check last files modified:\n";
system "ls -lt | head -6 ";

system "gmake depend" ;
system "gmake" ;

die "*** Compile failed!!! *** " unless ( -x $PROGRAM ) ;
open(MAKELOG,">Make.log");    # Over-write Make.log
print MAKELOG "$subv $NDX  $NDY  $dom_x0  $dom_y0  $dom_wx  $dom_wy \n" ;
close(MAKELOG);

if ( $COMPILE_ONLY) {     ## exit after make ##
    exit();
}


########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################

foreach $scenario ( @runs ) {
	print "STARTING RUN $scenario \n";

	$runlabel1    = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
	$runlabel2    = "${testv}_${scenario}_$year";   # NO SPACES! LONG (written into CDF files)

	#--- Change to RESDIR

	my $RESDIR = "$WORKDIR/$base/$scenario";
	system("mkdir -p $RESDIR");

	chdir "$RESDIR"; 

	# Erase any .nc files to avoid confusion:
	@nc_files = ( glob("$RESDIR/*.nc"),  glob("$RESDIR/*.nc.bz2") );
	foreach $f ( @nc_files) { unlink($f) };


#assign '-R';       # resets assignments    (huge)

@list_of_files = ();   # Keep list of data-files

for ($nnn = 1, $mm = $mm1; $mm <= $mm2; $mm++) {


    # Assign met files to fil001...etc.
    for ($n = 1; $n <= $month_days[$mm]; $n++) {
	# printf "file assign $mm1 $mm2 $mm %04d $n ==> ", $nnn;
	$nnn = &myfunc_mi2($n, $nnn, $mm);
	# printf "%04d $n\n", $nnn;
    }

    foreach $t ('snowc', 'natso2') {
	$old = sprintf "$SRCINPUTDIR/%s%02d.dat.170", $t, $mm;
	$new = sprintf "%s%02d.dat", $t, $mm;
        mylink( "Linking:", $old,$new ) ;
    }

    # Naming system simplified, 27/2/2003, ds
    #
    # Ozone,H2O2,OH, CH3COO2 data are needed for all model versions:

    foreach $bc ( qw ( D3_H2O2 D3_OH D3_CH3COO2 )) { #XD
            $old = sprintf "$BCDIR/$bc.%02d", $mm ;
            $new = sprintf "$bc.%02d", $mm ;
            mylink( "Linking $bc BCs:", $old,$new ) ;
    }

    foreach $bc ( qw ( D3_O3 )) { #XD
            $old = sprintf "$BCDIRO3/$bc.%02d", $mm ;
            $new = sprintf "$bc.%02d", $mm ;
            mylink( "Linking $bc BCs:", $old,$new ) ;
    }


#HF    $old = sprintf "$OZONEDIR/ozone.%02d", $mm ;
#HF    $new = sprintf "ozone.%02d", $mm ;
#HF    mylink( "Linking BCs:", $old,$new ) ;

#HF    $old = sprintf "$H2O2DIR/h2o2.%02d", $mm ;
#HF    $new = sprintf "h2o2.%02d", $mm ;
#HF    mylink( "Linking H2O2 BCs:", $old,$new ) ;

    $old = sprintf "$SRCINPUTDIR/lt21-nox.dat%02d", $mm;
    $new = sprintf "lightn%02d.dat", $mm;
    mylink( "Lightning : ", $old,$new ) ;
}


#-- individual treatment of the last record 

#BUG - FIX FOR 2000 NEEDED
if ( $NTERM > 100 ) {  # Cruide check that we aren't testing with NTERM=5
   if ( $mm2 == 12 ) {
        print "NEED TO SET H00 FROM NEXT YEAR \n";
	#97-fix $old = sprintf "$MetDir/f00.%04d0101", ($year+1)%100;
	$old = sprintf "$MetDir/f00.%04d0101", $year+1;
	$new = sprintf "fil%04d", $nnn;    #BUG_FIX : MISSING
        mylink( "LAST RECORD, NEED TO SET H00 FROM NEXT YEAR", $old,$new ) ;
   } else { #  Need 1st record of next month:
        $hhlast = 0 ;   #
        $ddlast = 1 ;
	$old = sprintf "$MetDir/f%02d.%04d%02d%02d", $hhlast, $year, $mm2+1, $ddlast;
	$new = sprintf "fil%04d", $nnn;
        mylink( "NEED TO SET H00 FROM NEXT MONTH", $old,$new ) ;
   }
}

#--- byteswap and assign the f* files

# double space for infield (io_infield=50, filxxx )
#assign '-F f77,bufa:96:2 p:fil%';

# double for outday, outmonth (io_out=80)
#assign '-F f77,bufa:96:2 p:out%';



#--- Prepare input data

# Emissions. This part is still a mixture, witk the old ko emission files
# left in for now. However, in future the only difference between
# MADE, MACHO and maade AERO-MADE should be in the numbre of pollutants
# used, which we can set in the @emisliost array.

# First, emission files are labelled e.g. gridSOx, whiuch we assign to
# emislist.sox to ensure compatability with the names (sox,...) used
# in the model.

%gridmap = ( "co" => "CO", "nh3" => "NH3", "voc" => "NMVOC", "sox" => "SOx",
             "nox" => "NOx" , "pm10" => "PM10", "pm25" => "PM25", "pmco" => "PMco" ) ;

# Emissions setup:

    #$old   = "$MyDataDir/femis.dat" ;
    $old   = "$femis_dir/femis.$scenario" ;
    die "ERROR::  NO femis $old !!!! \n" unless -r $old;
    $new   = "femis.dat";
    mylink( "Femis  ", $old,$new ) ;

    $old   = "$DaveDataDir/splits/pm25split.defaults.$Split" ;
    $new   = "pm25split.defaults";
    mylink( "Split pm25", $old,$new ) ;

    $old   = "$DaveDataDir/splits/vocsplit.defaults.$Split" ;
    $new   = "vocsplit.defaults";
    mylink( "Split voc", $old,$new ) ;

    $old   = "$DaveDataDir/splits/vocsplit.special.$Split" ;
    $new   = "vocsplit.special";
    mylink( "Split voc", $old,$new ) ;

foreach $poll  ( @emislist  ) {

   $old   = "$emisdir/grid$gridmap{$poll}" ;

   #
   #if ( ! $SR  && $poll =~ /pm/ ) {  # CRUDE FIX FOR NOW
   #    $old   = "$SVETLANA/Unify/MyData/emission/em2000/grid$gridmap{$poll}" ;
   #}
#
   $new   = "emislist.$poll";
   push @list_of_files, $new ;
   if( -r $new ) { unlink($new) };  # Get rid of any files from previous runs

# Africa and PM changes. 
# IMPORTANT. Do Africa first to get the NOx emissions needed for PM

$emissions_adjusted = $AFRICA_ADDED + $PM_ADDED;  # >=1 if something done

if ( $emissions_adjusted  ) {   # Copy since we want to change the file
   print "AFRICA $gridmap{$poll}\n";

   if ( $poll =~ /pm/ ) {
	print "Simply copy for  $new\n";
   	system("cat $old >  $new");  #ds use cat instead of copy to get write access
   } else {
	print "Add Africa for  $new\n";
	system("cat $old > tmp_$new");
	system("wc tmp_$new");
	system("cat tmp_$new $Africa/Africa$gridmap{$poll} > $new");
	system("wc $new");
	#system("cp $new africa_$new");  # For testing
	unlink("tmp_$new");
   }
} else {
   # Old system:
   print "No emissions adjustment\n";
   mylink( "Emis $poll : ", $old,$new ) ;
}

   $old   = "$timeseries/MonthlyFac.$poll" ;
   $new   = "MonthlyFac.$poll";
   mylink( "MonthlFac ", $old,$new ) ;

   $old   = "$timeseries/DailyFac.$poll" ;
   $new   = "DailyFac.$poll";
   mylink( "DailyFac ", $old,$new ) ;
}  # end of emissions poll loop

if ( $PM_ADDED ) {  # Add PM emissions based upon NOx inventory
   print "STARTING PM ADDITION\n";
   system("$DAVE/Unify/D_emis/mkp.pmemis_from_nox");
   #system("cat emislist.pm25 > test_emislist.pm25");
   #system("cat emislist.pmco > test_emislist.pmco");
}

# Surface measurement sites
# From $DAVE directory for these
    $old   = "$DataDir/sites.rep03" ;
    $new   =  "sites.dat";
    mylink("Sites ",  $old,$new ) ;

# Sondes
    $old   = "$DataDir/sondes.jan04" ;
    if( $SR          ) {$old   = "$DataDir/sondes.SR" };
    if( $MERLIN_CITY ) {$old   = "$DataDir/sondes.merlin_cities" };
    $new   = "sondes.dat";
    mylink( "Sondes", $old,$new ) ;

# Forest data
    #ds $old   = "$o3dir/forest.pcnt" ;
    #ds $new   = "forest.pcnt";
#MAR2004     $old   = "$DataDir/forests.tf2" ;
    $old   = "$DataDir/forests.mar2004b" ;
    $new   = "forest.dat";
    mylink( "Forest % cover", $old,$new ) ;

# Aircraft emissions
    $old   = "$o3dir/amilt42-nox.dat" ;
    $new   = "ancatmil.dat";
    mylink("Ancat ",  $old,$new ) ;

# Seasonal stuff
foreach $s ( keys(%seasons) ) {

    $old   = "$o3dir/a${s}t42-nox.dat" ;
    $new = sprintf "ancat$seasons{$s}.dat";
    mylink( "Ancat seasonal ", $old,$new ) ;

    $old   = "$o3dir/jclear.$s" ;
    $new = sprintf "jclear$seasons{$s}.dat";
    mylink( "Photolysis j-clear ", $old,$new ) ;

    $old   = "$o3dir/jcl1.$s" ;
    $new = sprintf "jcl1km$seasons{$s}.dat";
    mylink( "Photolysis j-1km ", $old,$new ) ;

    $old   = "$o3dir/jcl3.$s" ;
    $new = sprintf "jcl3km$seasons{$s}.dat";
    mylink( "Photolysis j-3km ", $old,$new ) ;
} 

    $old   = "$o3dir/rough.170" ;
    $new = sprintf "rough.170";
    mylink( "Roughness length", $old,$new ) ;

    #ds$old   = "$DataDir/landuse.nov2003" ;  #ds rv1_9_4 change
#MAR2004    $old   = "$DataDir/landuse.dec2003" ;  #ds rv1_9_4 change
    $old   = "$DataDir/landuse.mar2004" ;
    $new   = "landuse.dat";                #ds rv1_9_4 change
    mylink( "Landuse ", $old,$new ) ;

 # TMP LOCATION for some datafiles : MyDataDir
#ldefix foreach $datafile ( qw ( Volcanoes.dat tf2_gfac1.dat tf2_gfac2.dat tf2_biomass.dat ) ) {
#foreach $datafile ( qw ( Volcanoes.dat lde_gfac1.dat lde_gfac2.dat lde_biomass.dat ) ) {
foreach $datafile ( qw ( Volcanoes.dat  MM_gfac1.dat lde_gfac2.dat lde_biomass.dat ) ) {
    $old   = "$DataDir/$datafile" ;
    $new   = "$datafile" ;
    mylink( "$datafile", $old,$new ) ;
}


# FIX later - was the only emission control thingy....
  @exclus  = (9 ); #  NBOUND


#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------

 &flush(STDOUT);

# Link execultable also, since gridur is funny about these
# things....

  $LPROG = "prog.exe";
  mylink( "PROGRAM!!  ", $PROGRAM,$LPROG) ;

# Write out list of linked files to a shell-script, useful in case the program
# hangs or crashes:

 open(RMF,">Remove.sh");
 foreach $f ( @list_of_files ) {
     print RMF "rm $f \n";
     # print "REMOVE $f \n";
 }
 close(RMF);

# make file with input parameters (to be read by Unimod.f90)
system("rm INPUT.PARA");
open(TMP,">INPUT.PARA");
  print TMP "$NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n";
close(TMP);

foreach $exclu ( @exclus) {
    print "starting $PROGRAM with 
        NTERM $NTERM\nNASS $NASS\nEXCLU $exclu\nNDX $NDX\nNDY $NDY\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2";

# die #"STOP HERE  \n";

if ( $INTERACTIVE ) {
  die " -- INTERACTIVE: can now run for inputs: NTERM, NASS,  EXCLU, NDX, NDY 
  with mpirun -np 1 prog.exe << XXX
  $NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\nXXX"; 
}

    #open (PROG, "|mpirun -np $NPROC $PROGRAM") || 

    my $tmpdir = $ENV{PWD};
    $ENV{PWD} = $RESDIR;
    open (PROG, "|mpirun -np $NPROC $LPROG") || 
               die "Unable to execute $LPROG. Exiting.\\n" ;

    # open (PROG, "|mp
    #   open (PROG, "|$PROGRAM") || 
    #          die "Unable to execute $PROGRAM. Exiting.\\n";
    
    #print PROG "$NTERM\n\n$NASS\n$exclu\n$NDX\n$NDY\n\iyr_trend\n";
    #print PROG "$NTERM\n$NASS\n$exclu\n$NDX\n$NDY\n\iyr_trend\n";
#pw put into INPUT.PARA    print PROG "$NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n";
    close(PROG);
    $ENV{PWD} =$tmpdir;
}
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------

if ( -r core )  {
      #TEST unlink("core");   # Remove to save disk-space
      die "Error somewhere - Core dumped !!!!\n";
} else {
      #-- Done.
      print "\n  Eulmod: Successful exit at" . `date '+%Z %Y-%m-%d %T %j'` ." \n";
}
# Now clean up,
    foreach $f ( @list_of_files ) {
        unlink($f);
        #print "REMOVED $f \n";
    }
#move RunLog 
system("mv RunLog.out  ${runlabel1}_RunLog");

#tar sites and sondes. Use sondes to check as these are produced les frequently.
my $last_sondes = sprintf  "sondes.%02d%02d", $mm2, $yy;
print "LOOKING FOR LAST SITES $last_sondes\n";
if ( -r $last_sondes ) {
	print "FOUND LAST sondes $last_sondes\n";
	system("tar cvf ${runlabel1}.sites sites.*");
	system(	"tar cvf ${runlabel1}.sondes sondes.*");
}

# Make a list of big .nc files
   #@n_files = glob("$RESDIR/*_hour.nc");  # Takes too long, and only 5% saving
   @d_files = glob("$RESDIR/*_day.nc");
#    foreach $f ( @n_files, @d_files ) {
#        system("bzip2 -f $f");
#    }

#append to the list of files to be compressed
@bigfile_list = ( @bigfile_list, @n_files, @d_files );

################################## END OF SCENARIO RUNS ######################
}  ############################### END OF SCENARIO RUNS ######################
################################## END OF SCENARIO RUNS ######################

system "/usr/bin/ja -s"; #end of resource info gathering; give summary

exit 0;
# And compress the big files
$nproc_bsub = (split/\s+/,$ENV{'LSB_MCPU_HOSTS'})[1];

$Nnc=@bigfile_list;
$nleft = $Nnc % $nproc_bsub;

$nloops = ($Nnc-$nleft)/$nproc_bsub;
    print "number of files to be bunzip2:  $Nnc\n";
#exit;
#    print "number of processors: $nproc_bsub \n";
#    print "number of entire loops to be run:  $nloops\n";
#    print "number lefts: $nleft \n";

for ($i=1; $i<=$nloops; $i++){
    my $nkids=0;
    for ($j=1; $j<=$nproc_bsub; $j++){
	$nextfile=pop(@bigfile_list);
	
	if (my $pid = fork) {
	    $nkids++;
	    # parent, does nothing
#	    print "I am parent $pid \n";	    
	} else {
	    die "fork failed: $!" unless defined $pid;
	    # child, running compression
	    print "processing $nextfile \n";	    
	    system("bzip2", "-f", $nextfile);
	    exit 0; # child exits here
	}
    }
    if($nkids>0){
# parent wait:
    my $kid;
    do {
	# wait until all kids have finished
	$kid = waitpid(-1, WNOHANG);
    } until $kid > 0;
    print "We have just processed $i x $nproc_bsub files \n";	    
}
}

#prosess the files left
die "accounting error, $nleft $nproc_bsub" unless $nleft<$nproc_bsub;
    my $nkids=0;
    for ($j=1; $j<=$nleft; $j++){
	$nextfile=pop(@bigfile_list);
	
	if (my $pid = fork) {
	    $nkids++;
	    # parent, does nothing
	} else {
	    die "fork failed: $!" unless defined $pid;
	    # child, running compression
	    print "processing $nextfile \n";	    
	    system("bzip2", "-f", $nextfile);
	    exit 0; # child exits here
	}
    }
    if($nkids>0){
#parent is just waiting
    my $kid;
    do {
	# wait until all kids have finished
	$kid = waitpid(-1, WNOHANG);
    } until $kid > 0;
}



exit 0;


### SUBPROGRAMS ################################################################

sub leap_year {
    local ($y) = ($_[0]);

    if ($y < 20) {
	$y += 2000;
    } elsif ($y < 100) {
	$y += 1900;
    }

    if ($y % 400 == 0) {
	return 1;
    } elsif ($y % 100 == 0) {
	return 0;
    } else {
	return ($y % 4 == 0) ? 1 : 0;
    }
}

sub myfunc_mi2 {
    local ($dd,$nnn,$mm) = ($_[0], $_[1], $_[2]);
    local ($d,$hh,$old,$new);

    # meteorological data
    for ($hh = 0; $hh <= 21; $hh += 3) {
	$old = sprintf "$MetDir/f%02d.%04d%02d%02d", $hh, $year, $mm, $dd;
	$new = sprintf "fil%04d", $nnn;
        symlink $old,$new;
	push(@list_of_files , $new);    # For later deletion

	$nnn++;
    }

    return $nnn;
}


sub mylink {
  # links or assigns files from the original olcation (old) to
  # the new location (new) - generally the working directory.
  # Keeps track of all such linked files in list_of_files.
  # On gridur/odin we use symlink

    my ($text, $old,$new) = ($_[0], $_[1], $_[2]);

    symlink $old,$new || die "symlink $old $new failed : $!";
    # assign "-a $old   $new"  || die "assign $old $new failed : $!";

   print "$text $old => $new \n";
   push(@list_of_files , $new);    # For later deletion
}

sub calc_nterm {
  # Calculates the number of 3-hourly steps from month mm1 to
  #  mm2. Leap-years should already have been dealt with through
  #  the global month_days array which is used.
  #   $NTERM = 2921;    # works for non-leap year (365*8+1)

  local ($mm1,$mm2) = ($_[0], $_[1]) ;
  local($ndays=0,$i);

  foreach $i ( $mm1..$mm2 ) {
    $ndays += $month_days[$i] ;
  }
  $nterm = 1 + 8*$ndays ;
  print "Calculated NTERM = $nterm\n";
  return $nterm;
}

