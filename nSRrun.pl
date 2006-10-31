#!/usr/bin/perl

#Queue system commands start with # @ (these are not comments!)
#Be careful with uncommented lines starting with @ !

# @ job_name = emep
# @ job_type = parallel
# @ account_no = test001
## @ account_no       = nn2890k
# @ class = normal
# @ wall_clock_limit = 96:00:00
# @ node = 2
# @ tasks_per_node = 16
# @ node_usage = not_shared
# @ network.mpi = sn_all,shared,us
## @ total_tasks = 16
## @ max_processors =32
# @ resources        = ConsumableCpus(1) ConsumableMemory(300 mb)
## @ resources        = ConsumableCpus(16)
#
# @ notification     = never
# @ checkpoint       = no
# @ restart          = yes
#
# @ error            = job.$(Host).$(jobid).err
# @ output           = job.$(Host).$(jobid).out
# @ queue
#


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
# njord. The main advantage of this script is that
# the domain size  and input/output files and directories
# can be easily changed. The script does this by modifying the
# Par_ml.pat file to produce the corresponding
# Par_ml.f90 file with appropriate set-up.  
#
#    submit on njord from work directory with
#    >llsubmit frun.pl
#
#    check queue status with
#    >llq
#
#    kill with
#    >llcancel 3456
#
#
#
#
######################################################################

require "flush.pl";
my $MAKE = "make -f Makefile_njord";

#
@pattern_files = qw ( Par_ml.pat );


# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables 
#      are explained below, and derived variables set later.-

$year = "2000";
( $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

$SR = 1;   #Set to 1 for source-receptor stuff

$iyr_trend = $year;  
$iyr_trend = "2020" if $SR ;  # 2020 assumed for SR runs here

print "Year is $yy YEAR $year Trend year $iyr_trend\n";

$MetDir = "/work/emep/metdata/$year" ; 


#---  User-specific directories (changeable)

$PETER       = "mifapw/emep";      
$DAVE      = "mifads";      
$JOFFEN      = "mifajej";
$HILDE       = "mifahf";
$SVETLANA    = "mifast";
$HEIKO      = "mifahik";      
$ANNA      = "mifaab";      
$MICHAEL      = "mifaelg";      


$USER        =  $PETER ;      


my $HEMIS = 0;   #Set to 1 for Hemispheric run. Not possible yet
my $PM_ADDED     = 1;  # Adds in PM emissions from NOx inventory scaling
my $AFRICA_ADDED = 0;  # Adds in African emissions for y=1..11

$OZONE = "1"; 
$ACID = "0";     # Specify model type here
$SUM=$OZONE+$ACID;

  # check:
die "Must choose ACID or OZONE" if ( $OZONE+$ACID>1 or $OZONE+$ACID==0 );


# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

if ( $OZONE ) {
    @emislist = qw ( sox nox nh3 co voc pm25 pmco ); 
    $testv       = "rv2_6_1";
    # BC Logan for O3
    
} elsif ( $ACID ) {
    print "ACID not yet tested \n";	    
    
    exit;  # 
}

#Main general Data directory
$DataDir       = "/home/ntnu/mifapw/emep/Data";      

#Grid specific data
$DATA_LOCAL  = "/home/ntnu/mifapw/emep/Data/EMEP";

#User directories
$ProgDir     = "/home/ntnu/$USER/Unify/rv2_6_1";   # input of source-code
$WORKDIR     = "/work/mifapw/test";    # working and result directory
$MyDataDir   = "/home/ntnu/$USER/Unify/MyData";    # for each user's private input


#ds check: and change
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"  
#  unless $ENV{PWD} =~ /Unimod.$testv.$year/;


$Split       = "BASE_MAR2004" ;       
$NOxSplit       = "2000" ;               # Have CLE2020, MFR2020, 2000       
$Africa      = "$DATA_LOCAL/Africa";        # Emissions for Africa, y=1..11

$timeseries  = "$DataDir";

$version     = "Unimod" ;  
$PROGRAM     = "$ProgDir/$version";         # programme
$subv        = "$testv" ;                  # sub-version (to track changes)

#ds new for SR version. One lement arrays if SR not specified:
#$mk_SRfemis = "$DAVE/Unify/Progs/mkp.SRfemis";   #  Generates femis files

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later. 
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR
@runs        = ( $scenario );


#EMISSIONS
$EMIS_INP = "/home/ntnu/$HEIKO/Emissions/ModelInputModruns/Modrun06";
#$emisdir = "$EMIS_INP/2006-Trend2004-V7";
#$emisdir = "$EMIS_INP/2006_emis2020_SRlow_V7";
$emisdir = "$EMIS_INP/2006_emis2010_BL-E_V7";

$femis       = "$DataDir/femis.dat";      # emission control file


# Specify small domain if required. 
#                 x0   x1  y0   y1

@largedomain = (   1, 170,  1, 133 ) ;
##smalldomain = (  20, 167,  1, 122 ) ;    # OSPAR/HELCOM domain
##smalldomain = (  18, 169,  7, 124 ) ;     # OSPAR/HELCOM domain+border-south
@smalldomain = (  36, 167, 12, 122 ) ;    # EMEP domain
##@smalldomain = (  136, 167, 100, 122 ) ;      # (changeable)
##smalldomain = @largedomain ;     # If you want to run for the whole domain, 
		# simply uncomment this 


$RESET        = 0   ;  # usually 0 (false) is ok, but set to 1 for full restart
$COMPILE_ONLY = 0   ;  # usually 0 (false) is ok, but set to 1 for compile-only
$INTERACTIVE  = 0   ;  # usually 0 (false), but set to 1 to make program stop
# just before execution - so code can be run interactivel.

$NDX   = 8;           # Processors in x-direction
$NDY   = 4;           # Processors in y-direction
if ( $INTERACTIVE ) { $NDX = $NDY = 1 };


@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

$mm1   =  1;       # first month
$mm2   =  1;       # last month
$NTERM_CALC =  calc_nterm($mm1,$mm2);

$NTERM =   $NTERM_CALC;    # sets NTERM for whole time-period
# -- or --
#$NTERM = 6;       # for testing, simply reset here

print "NTERM_CALC = $NTERM_CALC, Used NTERM = $NTERM\n";

# <---------- end of normal use section ---------------------->
# <---------- start of SR   use section ---------------------->

#SR additions:
# Define all coutnries and nums here: Heiko's lib might be faster:

#SR additions:
# Define all coutnries and nums here: Heiko's lib might be faster:
%country_nums = (
  AL =>   1,  AT =>   2,  BE =>   3,  BG =>   4, FCS =>   5,
  DK =>   6,  FI =>   7,  FR =>   8,  # DE special: has number 60
  GR =>  11,  HU =>  12,  IS =>  13,  IE =>  14,  IT =>  15,
  LU =>  16,  NL =>  17,  NO =>  18,  PL =>  19,  PT =>  20,
  RO =>  21,  ES =>  22,  SE =>  23,  CH =>  24,  TR =>  25,
 FSU =>  26,  GB =>  27, REM =>  29, BAS =>  30,
 NOS =>  31, ATL =>  32, MED =>  33, BLS =>  34, NAT =>  35,
 RUO =>  36, RUP =>  37, RUA =>  38,  BY =>  39,  UA =>  40,
  MD =>  41, RUR =>  42,  EE =>  43,  LV =>  44,  LT =>  45,
  CZ =>  46,  SK =>  47,  SI =>  48,  HR =>  49,  BA =>  50,
  CS =>  51,  MK =>  52,  KZ =>  53,  GE =>  54,  CY =>  55,
  AM =>  56,  MT =>  57, ASI =>  58,  LI =>  59,  #special DE =>  60, RU =>  61,
  MC =>  62, NOA =>  63,  EU =>  64,  US =>  65,
  CA =>  66, BIC =>  67,  KG =>  68,  AZ =>  69,  DE =>  10, #DE special
  RU =>  42, #RU special
  RUX =>  71,  ATX =>  70
);

# EU countries:
@eu = qw ( AT BE DK FI FR DE GR IE IT NL PT ES SE GB LU );
@euaccess = qw ( HU PL CY CZ EE LT LV MT SK SI );
@eu25 = ( @eu, @euaccess );
@sea = qw ( NOS ATL MED BAS BLS );
@noneu = qw ( NO CH IS BG RO );
@emep = qw ( CS BY BA HR TR RU UA KZ MD MK GE AM AL AZ KG NOA ASI REM) ; 
@rains = qw ( CS BY BA HR TR RU UA MD MK AL NO CH BG RO ) ;  # RAINS only
@nonrains = qw ( KZ GE AM AZ KG NOA ASI REM );  # Not in RAINS, but in EMEP

@new =qw ( RUX   ATX ); 

@countries = qw ( IT ) ;  # List of wanted countries this run (at least 1)

@countries = @eu ;

@polls       = qw ( BASE );  # BASE, A, NP, V or S (at least 1)


$base        = "CLE2010";
$Split       = "CLE_MAR2004";    # VOC splits
$NOxSplit    = "CLE2020_ver2";    # NOx splits
$rednflag    = "P15";  # 15% reduction for label
$redn        = "0.85";  # 15% reduction



# <---------- end of user-changeable section ----------------->
#               (normally, that is...)



#--- Verify data directories
foreach $d (  $WORKDIR, $DATA_LOCAL, $DataDir,  $ProgDir) {
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

#--- Other eulmod configs

# 
# Check that we have an existing prog dir:
die "Wrong ProgDir: $ProgDir \n" unless -d $ProgDir;

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

#$nproc_bsub = (split/\s+/,$ENV{'LSB_MCPU_HOSTS'})[1];

#if ( ! $INTERACTIVE && ! $COMPILE_ONLY  && $NPROC != $nproc_bsub ) {
#    die " -- Requested wrong number of processors --
#      bsub asked for $nproc_bsub whereas NPROC = $NDX x $NDY = $NPROC \n";
#}




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

if ( $RESET ) { unlink ("Make.log") }  # Gone!
# ---- calculate domain widths
# (For the model, we need first x0, y0, then width (number of cells) in x and y)

$dom_x0 = $smalldomain[0]; 
$dom_y0 = $smalldomain[2]; 
$dom_wx = $smalldomain[1] - $smalldomain[0] + 1 ; 
$dom_wy = $smalldomain[3] - $smalldomain[2] + 1 ; 
$largdom_wx = $largedomain[1] - $largedomain[0] + 1 ; 
$largdom_wy = $largedomain[3] - $largedomain[2] + 1 ; 

open(MAKELOG,"<Make.log");
($oldversion, $olddx , $olddy , $old_x0, $old_y0, $old_wx, $old_wy, $old_lwx, $old_lwy) = split ' ', <MAKELOG> ;
close(MAKELOG);

print "From Make.log we had Subversion $subv
           Procs:  $olddx $olddy 
           Domain: $old_x0, $old_y0, $old_wx, $old_wy ,$old_lwx, $old_lwy \n";


if ( $oldversion ne $subv ) {
    
    print " We are changing version!!!!!............. 
            from $oldversion to $subv \n " ;
    $RESET = 1 ;
}

if ( $NDX      != $olddx  || $NDY      != $olddy  ||
     $dom_x0   != $old_x0 || $dom_y0   != $old_y0 ||
     $dom_wx   != $old_wx || $dom_wy   != $old_wy ||
     $largdom_wx   != $old_lwx || $largdom_wy   != $old_lwy       ) {
    
    print "Need to re-compile for new processor or domain  setup: 
           Procs:  $NDX $NDY
           Domain: $dom_x0, $dom_y0, $dom_wx, $dom_wy, $largdom_wx, $largdom_wy \n";
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
	    $line =~ s/largedomdx/$largdom_wx/ ;
	    $line =~ s/largedomdy/$largdom_wy/ ;
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

#to be sure that we don't use an old version (recommended while developing)
#unlink($PROGRAM);

system "$MAKE depend" ;
system "$MAKE" ;

#die "*** Compile failed!!! *** " unless ( -x $PROGRAM ) ;
open(MAKELOG,">Make.log");    # Over-write Make.log
print MAKELOG "$subv $NDX  $NDY  $dom_x0  $dom_y0  $dom_wx  $dom_wy $largdom_wx $largdom_wy \n" ;
close(MAKELOG);

if ( $COMPILE_ONLY) {     ## exit after make ##
    exit();
}




    
@list_of_files = ();   # Keep list of data-files



########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################


foreach my $pollut ( @polls ) {
	print "STARRTING POLLUT $pollut\n";
    foreach my $cc ( @countries ) {
	print "STARTING CC $cc\n";
    
    $ss = 100;#anthropogenic
    $scenario = "${base}_${cc}_${pollut}_${rednflag}" ;
    $scenario = "${base}" if $pollut eq "BASE" ;
    print "STARTING RUN $scenario \n";

    $runlabel1    = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
    $runlabel2    = "${testv}_${scenario}_$year";   # NO SPACES! LONG (written into CDF files)
    
    #--- Change to RESDIR
    
    my $RESDIR = "$WORKDIR/$year/$scenario";
    system("mkdir -p $RESDIR");
    
    chdir "$RESDIR"; 

#safety test:
#next if ( -r "$RESDIR/${runlabel1}_year.nc" );  
next if ( -r "$RESDIR/${runlabel1}_finished" );  

    
    
 
	my( $sox,$nox,$voc,$nh3,$testp,$co,$pm25,$pmco ) = ("1.0") x 8 ; # Initialise to 1.0
	if( $pollut eq "AV" ) { $voc = $nh3 = $redn  };
	if( $pollut eq "A" ) { $nh3 = $redn  };
	if( $pollut eq "V" ) { $voc = $redn  };
	if( $pollut eq "S" ) { $sox = $redn  };
	if( $pollut eq "N" ) { $nox = $redn  };
	if( $pollut eq "NP" ) { $nox = $pm25 = $pmco = $redn  };
	if( $pollut eq "SNP" ) { $sox = $nox = $pm25 = $pmco =  $redn  };
	if( $pollut eq "AN" ) { $nh3 = $nox =  $redn  };
	if( $pollut eq "SNAV" ) { $sox = $nox = $nh3 = $voc =  $redn  };
	#if( $pollut eq BASE ) then no change!

	my @ccloop = ( $cc );  		      # Normal, one-country option
##	@ccloop = @eu30  if( $cc eq "EU30" ); # Special multi-country option

	system("cat $DataDir/femis.dat > femis.dat");

	if ( $base eq "SN50" || $base eq "A50" || $base eq "M50" ){

	if($base eq "SN50") { 
	    if( $pollut eq "S" ) { $sox = $Xredn  };
	    if( $pollut eq "N" ) { $nox = $Xredn  };
	    if( $pollut eq "NP" ) { $nox = $Xredn  };
	    if( $pollut eq "SNP" ) { $sox = $nox = $Xredn  };
	    if( $pollut eq "AN" ) { $nox =  $Xredn  };
	    if( $pollut eq "SNAV" ) { $sox = $nox = $Xredn  };
	    }
	if($base eq "A50") { 
	    if( $pollut eq "AV" ) { $nh3 = $Xredn  };
	    if( $pollut eq "A" ) { $nh3 = $Xredn  };
	    if( $pollut eq "AN" ) { $nh3 =  $Xredn  };
	    if( $pollut eq "SNAV" ) { $nh3 = $Xredn  };
	}
	if($base eq "M50") { 
	    if( $pollut eq "AV" ) { $voc = $nh3 = $Xredn  };
	    if( $pollut eq "A" ) { $nh3 = $Xredn  };
	    if( $pollut eq "V" ) { $voc = $Xredn  };
	    if( $pollut eq "S" ) { $sox = $Xredn  };
	    if( $pollut eq "N" ) { $nox = $Xredn  };
	    if( $pollut eq "NP" ) { $nox = $pm25 = $pmco = $Xredn  };
	    if( $pollut eq "SNP" ) { $sox = $nox = $pm25 = $pmco =  $Xredn  };
	    if( $pollut eq "AN" ) { $nh3 = $nox =  $Xredn  };
	    if( $pollut eq "SNAV" ) { $sox = $nox = $nh3 = $voc =  $Xredn  };
	}

	foreach my $c ( @antropogenic ) {
#	my( $Xsox,$Xnox,$Xvoc,$Xnh3,$Xtestp,$Xco,$Xpm25,$Xpmco ) = ("0.5") x 8 ; # Initi
	my( $Xsox,$Xnox,$Xvoc,$Xnh3,$Xtestp,$Xco,$Xpm25,$Xpmco ) = ("1.0") x 8 ; # Initi
	$Xtestp= "1.0";
	if($base eq "SN50") { 
	    $Xsox = $Xnox = "0.5" ;
	    }
	if($base eq "A50") { 
	    $Xnh3 = "0.5" ;
	}
	if($base eq "M50") { 
	( $Xsox,$Xnox,$Xvoc,$Xnh3,$Xtestp,$Xco,$Xpm25,$Xpmco ) = ("0.5") x 8 ; # 
	}

	my $ccnum = $country_nums{$c} || die "ERROR!! No country Num!\n";
		system("echo \"$ccnum $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
		if ( $c eq "DE" ) {  # Add E. Germany too!
			system("echo \" 9  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
		    }
		if ( $c eq "ATL" ) {  # Add ATL outside EMEP
			system("echo \" 70  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
		    }
		if ( $c eq "RU" ) { # Add the other parts of Russia too! NB RU set to 42 (RUR)
			system("echo \" 38  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
			system("echo \" 37  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
			system("echo \" 36  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
			system("echo \" 71  $ss  $Xsox $Xnox $Xvoc $Xnh3 $Xtestp $Xco $Xpm25 $Xpmco\" >> femis.dat");
		}

        } # @antropogenic loop
    } #if X50

	foreach my $c ( @ccloop ) {
		my $ccnum = $country_nums{$c} || die "ERROR!! No country Num!\n";

		# NB still use sector 0 for AV runs. Should use 100 for others
		system("echo \"$ccnum $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
		if ( $c eq "DE" ) {  # Add E. Germany too!
			system("echo \" 9  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
		}
		if ( $c eq "ATL" ) {  # Add ATL outside EMEP
			system("echo \" 70  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
		}
		if ( $c eq "RU" ) { # Add the other parts of Russia too! NB RU set to 42 (RUR)
			system("echo \" 38  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
			system("echo \" 37  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
			system("echo \" 36  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
			system("echo \" 71  $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\" >> femis.dat");
		}



	} # ccloop
	system("echo NEW FEMIS:");
	system("cat femis.dat");


for ($nnn = 1, $mm = $mm1; $mm <= $mm2; $mm++) {
    
    # Assign met files to fil001...etc.
    for ($n = 1; $n <= $month_days[$mm]; $n++) {
	# printf "file assign $mm1 $mm2 $mm %04d $n ==> ", $nnn;
	$nnn = &myfunc_mi2($n, $nnn, $mm);
	# printf "%04d $n\n", $nnn;
	if($nnn > $NTERM){
	    last;
	}
    }
    
    
} #for ($nnn = 1


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
} #NTERM

    
    for ($nnn = 1, $mm = $mm1; $mm <= $mm2; $mm++) {
	

	foreach $t ('snowc', 'natso2') {
	    $old = sprintf "$DATA_LOCAL/%s%02d.dat.170", $t, $mm;
	    $new = sprintf "%s%02d.dat", $t, $mm;
	    mylink( "Linking:", $old,$new ) ;
	}
	
	
	$old = sprintf "$DataDir/lt21-nox.dat%02d", $mm;
	$new = sprintf "lightn%02d.dat", $mm;
	mylink( "Lightning : ", $old,$new ) ;
	
    } #for ($nnn = 1
    

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

#    $old   = "$DataDir/femis.dat" ;
#    $new   = "femis.dat";
#    mylink( "Femis  ", $old,$new ) ;

    $old   = "$DataDir/pm25split.defaults.$Split" ;
    $new   = "pm25split.defaults";
    mylink( "Split pm25", $old,$new ) ;

    $old   = "$DataDir/vocsplit.defaults.$Split" ;
    $new   = "vocsplit.defaults";
    mylink( "Split voc", $old,$new ) ;

    $old   = "$DataDir/vocsplit.special.$Split" ;
    $new   = "vocsplit.special";
    mylink( "Split voc", $old,$new ) ;

  #rv2_4_6  New, 30/5/2006, from Joffen

    $old   = "$DataDir/noxsplit.default.$NOxSplit" ;
    $new   = "noxsplit.defaults";
    mylink( "Split nox", $old,$new ) ;

    $old   = "$DataDir/noxsplit.special.$NOxSplit" ;
    $new   = "noxsplit.special";
    mylink( "Split nox", $old,$new ) ;


    foreach $poll  ( @emislist  ) {
	
	$old   = "$emisdir/grid$gridmap{$poll}" ;
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
	system("$DATA_LOCAL/emissions/mkp.pmemis_from_nox_ASI_NOA");
#	system("$DAVE/Unify/D_emis/mkp.pmemis_from_nox");
	#system("cat emislist.pm25 > test_emislist.pm25");
	#system("cat emislist.pmco > test_emislist.pmco");
    }
    
#BIC in NetCDF
    
    $old   = "$DATA_LOCAL/Boundary_and_Initial_Conditions.nc" ;
    $new   =  "Boundary_and_Initial_Conditions.nc";
    mylink("BIC ",  $old,$new ) ;
    
	
# Sondes
    $old   = "$DATA_LOCAL/sondes.dat" ;
    if( $SR ) {$old   = "$DATA_LOCAL/sondes.SR" };
    $new   = "sondes.dat";
    mylink( "Sondes", $old,$new ) ;

# Surface measurement sites
    $old   = "$DATA_LOCAL/sites.dat" ;
    $new   =  "sites.dat";
    mylink("Sites ",  $old,$new ) ;

    $old   = "$DATA_LOCAL/landuse.JUN06" ;
    ##$old   = "$DataDir/landuse.JUN06" ;
    $new   = "landuse.JUN06";
    mylink( "Landuse ", $old,$new ) ;

# Forest data
    $old   = "$DATA_LOCAL/forests.mar2004b" ;
    $new   = "forest.dat";
    mylink( "Forest % cover", $old,$new ) ;

# Aircraft emissions
    $old   = "$DataDir/amilt42-nox.dat" ;
    $new   = "ancatmil.dat";
    mylink("Ancat ",  $old,$new ) ;
    
# Seasonal stuff
    foreach $s ( keys(%seasons) ) {
	
	$old   = "$DataDir/a${s}t42-nox.dat" ;
	$new = sprintf "ancat$seasons{$s}.dat";
	mylink( "Ancat seasonal ", $old,$new ) ;
	
	$old   = "$DataDir/jclear.$s" ;
	$new = sprintf "jclear$seasons{$s}.dat";
	mylink( "Photolysis j-clear ", $old,$new ) ;
	
	$old   = "$DataDir/jcl1.$s" ;
	$new = sprintf "jcl1km$seasons{$s}.dat";
	mylink( "Photolysis j-1km ", $old,$new ) ;
	
	$old   = "$DataDir/jcl3.$s" ;
	$new = sprintf "jcl3km$seasons{$s}.dat";
	mylink( "Photolysis j-3km ", $old,$new ) ;
    } 
    
    $old   = "$DATA_LOCAL/rough.170" ;
    $new = sprintf "rough.170";
    mylink( "Roughness length", $old,$new ) ;
    
    #ds$old   = "$DataDir/landuse.nov2003" ;  #ds rv1_9_4 change
#MAR2004    $old   = "$DataDir/landuse.dec2003" ;  #ds rv1_9_4 change
    $old   = "$DATA_LOCAL/landuse.mar2004" ;
    $new   = "landuse.dat";                #ds rv1_9_4 change
    mylink( "Landuse ", $old,$new ) ;
    
    # TMP LOCATION for some datafiles : MyDataDir
    
    foreach $datafile ( qw ( Volcanoes.dat   ) ) {
	$old   = "$DATA_LOCAL/$datafile" ;
	$new   = "$datafile" ;
	mylink( "$datafile", $old,$new ) ;
    }

    foreach $datafile ( qw ( JUN06_gfac1.dat JUN06_gfac2.dat JUN06_biomass.dat ) ) { #JUN06
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
    
    ${startyear}=$year;
    ${startmonth}=$mm1;
    ${startday}=1;

$NASS  =  0;        # Set to one if "dump" of all concentrations wanted at end
    
# make file with input parameters (to be read by Unimod.f90)
    unlink("INPUT.PARA");
    open(TMP,">INPUT.PARA");
    print TMP "$NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n${startyear}\n${startmonth}\n${startday}\n";
    close(TMP);
    
    foreach $exclu ( @exclus) {
	print "starting $PROGRAM with 
        NTERM $NTERM\nNASS $NASS\nEXCLU $exclu\nNDX $NDX\nNDY $NDY\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2\n startyear ${startyear}\nstartmonth ${startmonth}\nstartday ${startday}\n";
	
# die #"STOP HERE  \n";
#    exit;
	
	
#if ( $INTERACTIVE ) {
#  die " -- INTERACTIVE: can now run for inputs: NTERM, NASS,  EXCLU, NDX, NDY 
#  with mpirun -np 1 prog.exe << XXX
#  $NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n${startyear}\n${startmonth}\n${startday}\nXXX"; 
#}


#	my $tmpdir = $ENV{PWD};
#	$ENV{PWD} = $RESDIR;
	

	system("pwd");
#    system("ls -l");
	


#system("ulimit -a");

	    open (PROG, "| ./$LPROG") || 
	    die "Unable to execute $LPROG. Exiting.\\n" ;
	close(PROG);
	
	
	
#	$ENV{PWD} =$tmpdir;
    } #foreach $exclu
            system("pwd");
#	exit;
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------

#system("rm -r core*");
    if ( -r core )  {
	#TEST unlink("core");   # Remove to save disk-space
	die "Error somewhere - Core dumped !!!!\n";
    } else {
	#-- Done.
	print "\n  Eulmod: Successful exit at" . `date '+%Z %Y-%m-%d %T %j'` ." \n";
    }
#move RunLog 
    system("mv RunLog.out  ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    system("echo Emis: $emisdir                 >> ${runlabel1}_RunLog");
    system("echo Version: $testv                >> ${runlabel1}_RunLog");
    system("echo Domain x0 $dom_x0 y0 $dom_y0 wx $dom_wx wy $dom_wy  >> ${runlabel1}_RunLog");
    system("echo Processors $NDX $NDY           >> ${runlabel1}_RunLog");
    system("echo Added? PM $PM_ADDED  Africa $AFRICA_ADDED >> ${runlabel1}_RunLog");
    system("echo SR?  $SR                       >> ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    system("echo femis: femis.$scenario         >> ${runlabel1}_RunLog");
    system("cat femis.dat >> ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    
    foreach $f ( @list_of_files ) {
        unlink($f);
        #print "REMOVED $f \n";
    }

    
#   exit;
    
#clean up work directories and links
    
    chdir "$RESDIR"; 
    foreach $f ( @list_of_files ) {
        unlink($f);
        #print "REMOVED $f \n";
    }
#    exit;
#tar sites and sondes. Use sondes to check as these are produced les frequently.
#    my $last_sondes = sprintf  "sondes.%02d%02d", $mm2, $yy;
#    print "LOOKING FOR LAST SITES $last_sondes\n";
#    if ( -r $last_sondes ) {
#	print "FOUND LAST sondes $last_sondes\n";
#	system("tar cvf ${runlabel1}.sites sites.*");
#	system(	"tar cvf ${runlabel1}.sondes sondes.*");
#    }
    

system(	"touch $RESDIR/${runlabel1}_finished");
last if $pollut eq "BASE"; 
    
################################## END OF SCENARIO RUNS ######################
}  # $cc ############################## END OF SCENARIO RUNS ######################
}  # $pollut ############################## END OF SCENARIO RUNS ######################
################################## END OF SCENARIO RUNS ######################


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

