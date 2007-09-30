#!/usr/bin/perl

#Common script for Njord and Snykov. Choose $NJORD=1 or $SNYKOV=1

#___________________________________________________________________
#Snykov queue commands

#Queue system commands start with #PBS (these are not comments!)
# lnodes= number of nodes, ppn=processor per node (max4)
#PBS -lnodes=6
# wall time limit of run 
#PBS -lwalltime=0:30:00
# lpmeme=memory to reserve per processor (max 4 or 16GB per node)
#PBS -lpmem=200MB
# account for billing
#PBS -A nn2890k
#___________________________________________________________________



#___________________________________________________________________
#Njord queue commands

#Queue system commands start with # @ (these are not comments!)
#Be careful with uncommented lines starting with @ !

# @ job_name = emep
# @ job_type = parallel
# @ account_no       = nn2890k
# @ class = express
# @ wall_clock_limit = 0:30:00
# @ node = 1
# @ tasks_per_node = 32
# @ node_usage = not_shared
# @ network.mpi = sn_all,shared,us
# @ resources        = ConsumableCpus(1) ConsumableMemory(300 mb)
# @ Requirements = ( Feature == "SMT" )
#
# @ notification     = never
# @ checkpoint       = no
# @ restart          = no
#
# @ error            = job.$(Host).$(jobid).err
# @ output           = job.$(Host).$(jobid).out
# @ queue
#
#___________________________________________________________________

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
#___________________________________________________________________
#
#    submit on SNOWSTORM from work directory with
#    >qsub srun.pl
#
#    check queue status with
#    >qstat -a
#
#    find more about job status (which nodes, when expected to start etc.)
#    >checkjob 3456
#
#    kill with
#    >qdel 3456
#
#
#
######################################################################

use 5.6.0;
use strict;
use warnings;
use File::Copy qw();

$| = 1; # autoflush STDOUT

#Choose one machine
my $NJORD=0; #1 if njord is used
my $SNYKOV=1; #1 if snykov (snowstorm) is used

# -j2 parallel make with 2 threads
my @MAKE = ("gmake",  "-j2", "--makefile=Makefile_snow");
  @MAKE = ("make", "-j2", "-f", "Makefile_njord") if $NJORD==1 ;
  die "Must choose SNYKOV **or** NJORD!\n" unless $NJORD+$SNYKOV==1;

my $SR= 0;     # Set to 1 if source-receptor calculation
              # check also variables in package EMEP::Sr below!!
 
# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables 
#      are explained below, and derived variables set later.-

my $year = "2002";
( my $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

my $iyr_trend = $year;
$iyr_trend = "2020" if $SR ;  # 2020 assumed for SR runs here

print "Year is $yy YEAR $year Trend year $iyr_trend\n";


#---  User-specific directories (changeable)

my $PETER      = "mifapw/emep";      
my $DAVE       = "mifads";      
my $JOFFEN     = "mifajej";
my $HILDE      = "mifahf";
my $SVETLANA   = "mifast";
my $HEIKO      = "mifahik";      
my $ANNA       = "mifaab";      
my $MICHAEL    = "michaelg";      
my $SEMEENA    = "mifasv";      


my $USER        =  $DAVE;
my ($HOMEROOT, $WORKROOT, $MetDir);
our $DataDir;
if ($SNYKOV){
    $HOMEROOT       = "/home";      
    $WORKROOT     = "/global/work";      
    $DataDir       = "/home/mifapw/emep_common/Data";
    $MetDir        = "$DataDir/EMEP/metdata/$year" ;
} else {
    $HOMEROOT       = "/home/ntnu";      
    $WORKROOT     = "/work";      
    $MetDir        = "/work/emep/metdata/$year" ;
    $DataDir       = "/home/ntnu/mifapw/emep/Data";
}
# DataDir    = Main general Data directory

my $DATA_LOCAL    = "$DataDir/EMEP";    # Grid specific data
my $NDATA_LOCAL   = "/home/mifads/Unify/Data/EMEP";  #  TMP! New Input style


my $HEMIS = 0;   #Set to 1 for Hemispheric run. Not possible yet
my $PM_ADDED     = 0;  # Adds in PM emissions from NOx inventory scaling
my $AFRICA_ADDED = 0;  # Adds in African emissions for y=1..11

my $OZONE = "1"; 
my $ACID = "0";     # Specify model type here, and check:
  die "Must choose ACID or OZONE" unless $OZONE+$ACID==1;


# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

my (@emislist, $testv);
if ( $OZONE ) {
    @emislist = qw ( sox nox nh3 co voc pm25 pmco ); 
    $testv       = "rv3beta30";
    
} elsif ( $ACID ) {
    die "ACID not yet tested \n";	    
}


#User directories
my $ProgDir     = "$HOMEROOT/$USER/Unify/Unimod.$testv";   # input of source-code
my $WORKDIR     = "$WORKROOT/$USER/$testv.$year";    # working and result directory
my $MyDataDir   = "$HOMEROOT/$USER/Unify/MyData";    # for each user's private input

#ds check: and change
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"  
#  unless $ENV{PWD} =~ /Unimod.$testv.$year/;


my $Split       = "BASE_MAR2004" ;       
my $NOxSplit       = "2000" ;               # Have CLE2020, MFR2020, 2000       
#$NOxSplit    = "CLE2020_ver2";    # NOx splits
my $Africa      = "$DATA_LOCAL/Africa";        # Emissions for Africa, y=1..11

my $timeseries  = "$DataDir";

my $version     = "Unimod" ;  
my $PROGRAM     = "$ProgDir/$version";         # programme
my $subv        = "$testv" ;                  # sub-version (to track changes)

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later. 
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR
my @runs        = ( $scenario );


#EMISSIONS
my $EMIS_INP = "$HOMEROOT/$HEIKO/Emissions/ModelInputModruns/Modrun06";
#my $emisdir = "$EMIS_INP/2006_emis2010_BL-E_V7";
my $emisdir = "$EMIS_INP/2006-Trend${year}-V7";
$emisdir = "$EMIS_INP/2006-Trend2004-V7" if $year > 2004;
my $pm_emisdir = $emisdir;
$pm_emisdir = "$EMIS_INP/2006-Trend2000-V7"  if $year < 2000;
 

my $RESET        = 0 ;  # usually 0 (false) is ok, but set to 1 for full restart
my $COMPILE_ONLY = 0 ;  # usually 0 (false) is ok, but set to 1 for compile-only
my $INTERACTIVE  = 0 ;  # usually 0 (false), but set to 1 to make program stop
my $DRY_RUN      = 0 ;  # Test script without running model (but compiling)

# just before execution - so code can be run interactivel.

# NDX, NDY  now set in ModelConstants_ml - use perl to extract these
# values and check against submission:

open(IN,"<$ProgDir/ModelConstants_ml.f90");
my ( $NDX, $NDY ); # Processors in x-, y-, direction
while(my $line = <IN>){
    $line=~ s/!.*//; # Get rid of comment lines
    $NDX = $1 if $line =~ /\W+ NPROCX \W+ (\d+) /x ;
    $NDY = $1 if $line =~ /\W+ NPROCY \W+ (\d+) /x ;
}
close(IN);
my $NPROC =  $NDX * $NDY ;
print "ModelConstants has: NDX = $NDX NDY = $NDY  =>  NPROC = $NPROC\n";

if ( $ENV{PBS_NODEFILE} ) {
   $_ =  `wc -l $ENV{PBS_NODEFILE}`;
   my $RUN_NPROC;
   ( $RUN_NPROC, undef ) = split;
   print "Qsub has: lnodes $RUN_NPROC\n";
   die "Error: Wrong number of lnodes!\n" unless $NPROC == $RUN_NPROC;
} else {
   print "skip nodefile check on interactive runs\n";
}

# P B S -lnodes=$NPROC
#(Didn't work: my $NCPUS=`wc -l $ENV{PBS_NODEFILE} | awk '{print $1}'`;)



my @month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

my $mm1   =  "06";       # first month, use 2-digits!
my $mm2   =  "06";       # last month, use 2-digits!
my $NTERM_CALC =  calc_nterm($mm1,$mm2);

my $NTERM =   $NTERM_CALC;    # sets NTERM for whole time-period
# -- or --
$NTERM = 2;       # for testing, simply reset here

print "NTERM_CALC = $NTERM_CALC, Used NTERM = $NTERM\n";

# <---------- end of normal use section ---------------------->
# <---------- end of user-changeable section ----------------->
#               (normally, that is...)

if ($SR) {
    print "SR is true\n";
    @runs = EMEP::Sr::initRuns();
}


#--- Verify data directories
mkdir_p($WORKDIR);
foreach my $d (  $WORKDIR, $DATA_LOCAL, $DataDir,  $ProgDir) {
    unless ( -d "$d" &&  -x _ && -r _ ) {
	die "*** ERROR *** directory $d not accessible. Exiting.\n";
    }
}


#--- Other eulmod configs

# 
# Check that we have an existing prog dir:
die "Wrong ProgDir: $ProgDir \n" unless -d $ProgDir;


# and check that we have the same from the bsub command, by
# accessing the LSB_MCPU_HOSTS environment variable.
# (For a 5 cpu job this looks like: LSB_MCPU_HOSTS=gridur.ntnu.no 5)

#$nproc_bsub = (split/\s+/,$ENV{'LSB_MCPU_HOSTS'})[1];

#if ( ! $INTERACTIVE && ! $COMPILE_ONLY  && $NPROC != $nproc_bsub ) {
#    die " -- Requested wrong number of processors --
#      bsub asked for $nproc_bsub whereas NPROC = $NDX x $NDY = $NPROC \n";
#}




#--- Calendar stuff
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);


#--- adjust for leap year
$month_days[2] += leap_year($year);

# --- Start new compilations if needed

# --- We check if the already-compiled version of $PROGRAM is the same
#     as the one we are asking for. The old configuration should have 
#     been stored in Make.log.  Read in this and split terms into 
#     : model, emis, ndx, ndy and domain sizes


chdir "$ProgDir";

#-- generate Makefile each time, to avoid forgetting changed "pat" file!


if ( $RESET ) { ########## Recompile everything!
    
    # For now, we simply recompile everything!
    system(@MAKE, "clean");
    system(@MAKE, "depend");
    system(@MAKE, "all");
}
system "pwd";
print "Check last files modified:\n";
system "ls -lt | head -6 ";

#to be sure that we don't use an old version (recommended while developing)
#unlink($PROGRAM);

system (@MAKE, "depend") ;
system (@MAKE, "all");

die "Done. COMPILE ONLY\n" if  $COMPILE_ONLY;  ## exit after make ##

    
my @list_of_files = ();   # Keep list of data-files



########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################


foreach my $scenflag ( @runs ) {
    if ($SR) {
	$scenario = EMEP::Sr::getScenario($scenflag);
    } else {
	$scenario = $scenflag;
    }
    print "STARTING RUN $scenario \n";

    my $runlabel1    = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
    my $runlabel2    = "${testv}_${scenario}_$year\_$iyr_trend";   # NO SPACES! LONG (written into CDF files)

    my $RESDIR = "$WORKDIR/$scenario";
    mkdir_p($RESDIR);

    chdir $RESDIR;   ############ ------ Change to RESDIR
    print "Working in directory: $RESDIR\n";

    my $METformat="felt";
    if ($METformat eq "felt") {
	my $nnn = 1;
	for (my $mm = $mm1; $mm <= $mm2; $mm++) {
	    
	    # Assign met files to fil001...etc.
	    for (my $n = 1; $n <= $month_days[$mm]; $n++) {
		$nnn = metlink($n, $nnn, $mm);
		last if $nnn > $NTERM;
	    }
	}
	
	my $mmlast = $mm2 + 1;
	my $yylast = $year;
	if ( $mmlast > 12 && $NTERM > 200 ) { # Crude check that we aren't testing with NTERM=5
	    $yylast = $yylast + 1;
	    $mmlast = 1;
	}
	my $old = sprintf "$MetDir/f00.%04d%02d01", $yylast, $mmlast;
	my $new = sprintf "fil%04d", $nnn;
	mylink( "LAST RECORD SET: ", $old,$new ) ;

    }else{

#meteorology in NetCDF	
	for (my $mm = $mm1; $mm <= $mm2; $mm++) {
	    for (my $n = 1; $n <= $month_days[$mm]; $n++) {
		my  $old = sprintf "$MetDir/meteo${year}%02d%02d.nc", $mm,$n;
		my  $new = sprintf "meteo${year}%02d%02d.nc", $mm,$n;
		mylink( "Linking:", $old,$new ) ;
	    }   
	}
	my $mmlast = $mm2 + 1;
	my $yylast = $year;
	if ( $mmlast > 12 && $NTERM > 200 ) { # Crude check that we aren't testing with NTERM=5
	    $yylast = $yylast + 1;
	    $mmlast = 1;
	}
	my   $old = sprintf "$MetDir/meteo%02d%02d01.nc", $yylast, $mmlast;
	my   $new = sprintf "meteo%02d%02d01.nc", $yylast, $mmlast;
	mylink( "LAST RECORD SET: ", $old,$new ) ;
    }

#=================== INPUT FILES =========================================
# ToDo Change noxsplit.default to defaults, as with voc (also in Unimod)
#   AFRICA, PM_ADDED --- fix elsewhere.

    my %ifile   = ();   # List of input data-files

# First, emission files are labelled e.g. gridSOx, whiuch we assign to
# emislist.sox to ensure compatability with the names (sox,...) used
# in the model.

my %gridmap = ( "co" => "CO", "nh3" => "NH3", "voc" => "NMVOC", "sox" => "SOx",
	 "nox" => "NOx" , "pm10" => "PM10", "pm25" => "PM25", "pmco" => "PMco" ) ;

    foreach my $poll  ( @emislist  ) {
        my $dir = $emisdir;
        $dir = $pm_emisdir if $poll =~ /pm/;   # FIX needed prior to 2000
        $ifile{"$dir/grid$gridmap{$poll}"} = "emislist.$poll";
        $ifile{"$timeseries/MonthlyFac.$poll"} = "MonthlyFac.$poll";
        $ifile{"$timeseries/DailyFac.$poll"} = "DailyFac.$poll";
    }
    
    foreach my $mmm ( $mm1  .. $mm2 ) {
	my $mm = sprintf "%2.2d", $mmm ; # WHY DO WE NEED THIS?????
	$ifile{"$DATA_LOCAL/snowc$mm.dat.170"} =  "snowc$mm.dat";
	$ifile{"$DATA_LOCAL/natso2$mm.dat.170"} =  "natso2$mm.dat";
	$ifile{"$DataDir/lt21-nox.dat$mm"} =  "lightn$mm.dat";
    }

# Emissions setup:
    $ifile{"$DataDir/femis.dat"} =  "femis.dat";
    $ifile{"$DataDir/pm25split.defaults.$Split"} = "pm25split.defaults";
    $ifile{"$DataDir/vocsplit.defaults.$Split"} = "vocsplit.defaults";
    $ifile{"$DataDir/vocsplit.special.$Split"} = "vocsplit.special";
    $ifile{"$DataDir/noxsplit.default.$NOxSplit"} = "noxsplit.defaults"; # NB - no
    $ifile{"$DataDir/noxsplit.special.$Split"} = "noxsplit.special";
    $ifile{"$DATA_LOCAL/Boundary_and_Initial_Conditions.nc"} =
                                   "Boundary_and_Initial_Conditions.nc";
    $ifile{"$DataDir/amilt42-nox.dat"} = "ancatmil.dat";#RENAME TO AIRCARAFT?!

  # new inputs style (Aug 2007)  with compulsory headers:
    $ifile{"$NDATA_LOCAL/Inputs.BVOC"} = "Inputs.BVOC";
    $ifile{"$NDATA_LOCAL/Inputs.Landuse"} = "Inputs.Landuse";
    $ifile{"$NDATA_LOCAL/Inputs.LandDefs"} = "Inputs.LandDefs";
    $ifile{"$NDATA_LOCAL/Inputs_DO3SE.csv"} = "Inputs_DO3SE.csv";
    $ifile{"$NDATA_LOCAL/sites.dat"} = "sites.dat";
    $ifile{"$NDATA_LOCAL/sondes.dat"} = "sondes.dat";
  
    
# Seasonal stuff  ----    Can't we improve this? e.g. every month?
    my %seasons = ( "jan" => "01", "apr" => "02", "jul" => "03" , "oct"=> "04");

    foreach my $s ( keys(%seasons) ) {
	$ifile{"$DataDir/a${s}t42-nox.dat"} = "ancat$seasons{$s}.dat";
	$ifile{"$DataDir/jclear.$s"} = "jclear$seasons{$s}.dat";
	$ifile{"$DataDir/jcl1.$s"} = "jcl1km$seasons{$s}.dat";
	$ifile{"$DataDir/jcl3.$s"} = "jcl3km$seasons{$s}.dat";
    } 
    
    $ifile{"$DATA_LOCAL/rough.170"} = "rough.170"; # Roughness length;
    $ifile{"$DATA_LOCAL/Volcanoes.dat"} = "Volcanoes.dat";


    foreach my $old ( sort keys %ifile ) {  # CHECK and LINK
	if ( -r $old ) {
		my $new =  $ifile{$old};
       		mylink( "Inputs: ", $old,$new ) ;
	} else {
        	print "Missing Input $old !!!\n";
		die "ERROR: Missing OLD $old\n" unless $old =~ /special/;
	}
    }

    if ($SR) {
	EMEP::Sr::generate_updated_femis(@$scenflag);
    }

#=================== INPUT FILES =========================================

    
# FIX later - was the only emission control thingy....
    my @exclus  = (9 ); #  NBOUND
    
    
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
    
    print "\n";
    
# Link execultable also, since gridur is funny about these
# things....
    
    my $LPROG = "prog.exe";
    #mylink( "PROGRAM!!  ", $PROGRAM,$LPROG) ;
    cp ($PROGRAM, $LPROG) or die "cannot copy $PROGRAM to $LPROG\n";
    push(@list_of_files , $LPROG);    # For later deletion
    
# Write out list of linked files to a shell-script, useful in case the program
# hangs or crashes:
    
    open(RMF,">Remove.sh");
    foreach my $f ( @list_of_files ) {	print RMF "rm $f \n"; }
    print RMF "rm $LPROG\n";   # Delete executable also
    close(RMF);
    
    my $startyear = $year;
    my $startmonth = $mm1;
    my $startday = 1;

    my $NASS  =  0;   # Set to one if "dump" of all concentrations wanted at end
    
# make file with input parameters (to be read by Unimod.f90)
    unlink("INPUT.PARA");
    open(TMP,">INPUT.PARA");
    print TMP "$NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n${startyear}\n${startmonth}\n${startday}\n";
    close(TMP);
    
    foreach my $exclu ( @exclus) {
	print "starting $PROGRAM with 
        NTERM $NTERM\nNASS $NASS\nEXCLU $exclu\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2\n startyear ${startyear}\nstartmonth ${startmonth}\nstartday ${startday}\n";
	

	my $PRERUN = "";
	$PRERUN = "scampiexec " if $SNYKOV;
	if ($DRY_RUN) {
	    print "DRY_RUN: not running '| $PRERUN ./$LPROG'\n";
	} else {
	    open (PROG, "| $PRERUN ./$LPROG") || 
		die "Unable to execute $LPROG. Exiting.\\n" ;
	    close(PROG);
	}
	
    } #foreach $exclu
            system("pwd");
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------
#------------    End of Run model -------------------------------------

    if ( -r "core" )  {
	die "Error somewhere - Core dumped !!!!\n";
    } elsif ( -r "Timing.out" ) { #-- Done  :-)
        print "\n  Eulmod: Successful exit at" . `date '+%Z %Y-%m-%d %T %j'` ." \n";
    } else {
        print "\n  The program stopped abnormally!! \n" unless $DRY_RUN;
    }

#move RunLog 
    rename "RunLog.out",  "${runlabel1}_RunLog"
	or warn "cannot mv RunLog.out ${runlabel1}_RunLog\n" unless $DRY_RUN;
    open RUNLOG, ">> ${runlabel1}_RunLog"
	or die "cannot append ${runlabel1}_RunLog: $!\n";
    print RUNLOG <<"EOT";
Emission units: Gg/year       
------------------------------
Emissions: $emisdir           
Version: $testv               
Processors $NDX $NDY          	
Added? PM $PM_ADDED  Africa $AFRICA_ADDED
SR?  $SR                      
iyr_trend: $iyr_trend         
------------------------------
femis: femis.$scenario        
------------------------------
EOT
    close RUNLOG;

#clean up work directories and links   
    if ($DRY_RUN) { # keep femis.dat
	@list_of_files = grep {$_ ne 'femis.dat'} @list_of_files;
    }
    unlink ( @list_of_files );

#tar sites and sondes. Use sondes to check as these are produced les frequently.
    my $last_sondes = sprintf  "sondes.%02d%02d", $mm2, $yy;
    print "LOOKING FOR LAST SITES $last_sondes\n";
    if ( -r $last_sondes ) {
	print "FOUND LAST sondes $last_sondes\n";
	system("tar cvf ${runlabel1}.sites sites.*");
	system("tar cvf ${runlabel1}.sondes sondes.*");
    }


################################## END OF RUNS ######################
}  ############################### END OF RUNS ######################
################################## END OF RUNS ######################


exit 0;


### SUBPROGRAMS ################################################################

sub leap_year {
    my ($y) = ($_[0]);

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

sub metlink {  #---- meteorological data
    my ($dd,$nnn,$mm) = ($_[0], $_[1], $_[2]);

    for (my $hh = 0; $hh <= 21; $hh += 3) {
	my $old = sprintf "$MetDir/f%02d.%04d%02d%02d", $hh, $year, $mm, $dd;
	my $new = sprintf "fil%04d", $nnn;
	mylink("Met:", $old, $new);
	$nnn++;
    }
    return $nnn;
}


sub mylink {
  # links files from the original olcation (old) to
  # the new location (new) - generally the working directory.
  # Keeps track of all such linked files in list_of_files.

    my ($text, $old,$new) = ($_[0], $_[1], $_[2]);

    symlink $old,$new || die "symlink $old $new failed : $!";

   print "$text $old => $new \n";
   push(@list_of_files , $new);    # For later deletion
}

sub touch {
    # simple touch -c implementation
    my (@fileGlobs) = @_;
    my @files;
    foreach my $fileGlob (@fileGlobs) {
	push @files, glob($fileGlob);
    }
    utime undef, undef, @files;
}

sub cp {
    # copy, preserving permissions (stupid File::Copy::cp does not)
    my ($from, $to, @extraArgs) = @_;
    my $retVal = File::Copy::cp($from, $to, @extraArgs);
    my $perm = (stat $from)[2] & 07777;
    chmod($perm, $to);
    return $retVal;
}

sub mkdir_p {
    # mkdir -p on unix platforms
    # does NOT fail on existing directories!
    my ($dir) = @_;
    $dir =~ s:/$::; # remove final /
    my $curdir = './';
    if ($dir =~ s:^/::) {
	$curdir = '/';
    }
    my @path = split ('/', $dir);
    while (my $next = shift(@path)) {
	$curdir .= $next . '/';
	if (! -d $curdir) {
	    mkdir $curdir or die "cannot mkdir $curdir: $!\n";
	}
    }
    return 1;
}

sub calc_nterm {
  # Calculates the number of 3-hourly steps from month mm1 to
  #  mm2. Leap-years should already have been dealt with through
  #  the global month_days array which is used.
  #   $NTERM = 2921;    # works for non-leap year (365*8+1)

  my ($mm1,$mm2) = ($_[0], $_[1]) ;
  my $ndays=0;
  foreach my $i ( $mm1..$mm2 ) {
    $ndays += $month_days[$i] ;
  }
  my $nterm = 1 + 8*$ndays ;

  print "Calculated NTERM = $nterm\n";
  return $nterm;
}


##############################################################
##### Stuff for Source receptor matrisses             ########
##############################################################
package EMEP::Sr;

my (%country_nums, @eu15, @euNew04, @eu25, @euNew06, @eu27, @sea, @noneu, @emep, @external);
our ($base, $Split, $NOxSplit, $rednflag, $redn, @countries, @polls);

INIT {
########################################
# Define all countries and nums here: ##
########################################
%country_nums = (
  AL =>   1,  AT =>   2,  BE =>   3,  BG =>   4, FCS =>   5,
  DK =>   6,  FI =>   7,  FR =>   8,
  GR =>  11,  HU =>  12,  IS =>  13,  IE =>  14,  IT =>  15,
  LU =>  16,  NL =>  17,  NO =>  18,  PL =>  19,  PT =>  20,
  RO =>  21,  ES =>  22,  SE =>  23,  CH =>  24,  TR =>  25,
 FSU =>  26,  GB =>  27, REM =>  29, BAS =>  30,
 NOS =>  31, ATL =>  32, MED =>  33, BLS =>  34, NAT =>  35,
 RUO =>  36, RUP =>  37, RUA =>  38,  BY =>  39,  UA =>  40,
  MD =>  41, RUR =>  42,  EE =>  43,  LV =>  44,  LT =>  45,
  CZ =>  46,  SK =>  47,  SI =>  48,  HR =>  49,  BA =>  50,
  RS => 72,   ME =>  73,  MK =>  52,  KZ =>  53,  GE =>  54,  CY =>  55,
  AM =>  56,  MT =>  57, ASI =>  58,  LI =>  59,  DE =>  60, RU =>  61,
  MC =>  62, NOA =>  63,  EU =>  64,  US =>  65,
  CA =>  66, BIC =>  67,  KG =>  68,  AZ =>  69,
  RUX =>  71,  ATX =>  70,
 BA2 => 302, BA3 => 303, BA4 => 304, BA5 => 305, BA6 => 306, # Baltic sep.
 BA7 => 307, BA8 => 308, BA9 => 309,
 NS2 => 312, NS3 => 313, NS4 => 314, NS5 => 315, NS6 => 316, # N. Sea sep.
 NS7 => 317, NS8 => 318, NS9 => 319,
 AT2 => 322, AT3 => 323, AT4 => 324, AT5 => 325, AT6 => 326, # Atlant. sep.
 AT7 => 327, AT8 => 328, AT9 => 329,
 ME2 => 332, ME3 => 333, ME4 => 334, ME5 => 335, ME6 => 336, # Medit. sep.
 ME7 => 337, ME8 => 338, ME9 => 339,
 BL2 => 342, BL3 => 343, BL4 => 344, BL5 => 345, BL6 => 346, # Bl. Sea sep.
 BL7 => 347, BL8 => 348, BL9 => 349
);

# EU countries:
@eu15 = qw ( AT BE DK FI FR DE GR IE IT NL PT ES SE GB LU );
@euNew04 = qw ( HU PL CY CZ EE LT LV MT SK SI );
@eu25 = ( @eu15, @euNew04 );
@euNew06 = qw(BG RO);
@eu27 = (@eu25, @euNew06);
@sea = qw ( NOS ATL MED BAS BLS );
@noneu = qw ( NO CH IS );
@emep = qw ( RS ME BY BA HR TR RU UA KZ MD MK GE AM AL AZ KG NOA ASI REM) ; 
@external =qw ( RUX   ATX ); 

########################################
# End of country definitions          ##
########################################


################################
#### start of SR parameters ####
################################
$base        = "CLE";
$Split       = "CLE_MAR2004";     # Defualt (IER-based) VOC splits
$NOxSplit    = "CLE2020_ver2";    # Default scenario (IER-based) VOC splits
$rednflag    = "P15";  # 15% reduction label
$redn        = "0.85"; # 15% reduction

# modify those to fill up your queues for SR effectively!!!
@countries  = (@eu27, @sea, @noneu, @emep);
@polls       = qw (BASE NP A V S );  #  (any, all, at least 1)
################################
#### end of SR parameters   ####
################################
}


sub initRuns {
    my @runs;
    foreach my $cc (@countries) {
	foreach my $poll (@polls) {
	    push @runs, [$cc, $poll, $redn];
	    if ($poll eq 'BASE') {
		# run BASE only once (for exactly one cc)!!!
		@polls = grep {'BASE' ne $_} @polls;
	    }
	}
    }
    return @runs;
}


sub getScenario {
    my ($scenflag) = @_;
    my $cc = $scenflag->[0];
    my $pollut = $scenflag->[1];
    my $scenario = "${base}_${cc}_${pollut}_${rednflag}";
    $scenario = "${base}" if $pollut eq "BASE" ;
    return $scenario;
}

sub generate_updated_femis {
    my ($cc, $pollut, $redn) = @_;
   # Initialise to 1.0:
    my( $sox,$nox,$voc,$nh3,$testp,$co,$pm25,$pmco ) = ("1.0") x 8 ; 
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

    my $femisdat = slurp("$DataDir/femis.dat");

    my $ccnum = $country_nums{$cc} || die "ERROR!! No country Num for $cc!\n";

    # using 0 here as long as emissions are guaranteed to contain either
    # only anthropogenic or only natural emissions perl 'country'
    my $ss = 0; # 100 = antropogenic sectors (1..10)
                # 0 = all sectors
    $femisdat .= "$ccnum $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    if ( $cc eq "DE" ) {  # Add splitted countries
	foreach my $cx (9, 10) {
	    $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
	}
    }
    if ( $cc eq "RU" ) { # Add splitted and external RU
	foreach my $cx (36..38, 42, 71) {
	    $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
	}
    }
    if ( $cc eq "ATL" ) {  # Add ATL outside EMEP
	foreach my $cx (70) {
	    $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
	}
    }
    if ( 30 <= $ccnum and $ccnum <= 34) { # add splitted sea areas
	for (my $cx = 10 * $ccnum + 2; $cx <= 10 * $ccnum + 9; $cx++) {
	    $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
	}
    }
    unlink "femis.dat" if -l "femis.dat";
    open FEMIS, ">femis.dat" or die "Cannot write femis.dat: $!\n";
    print FEMIS $femisdat;
    close FEMIS;

    # and to the logfile
    print "NEW FEMIS\n", $femisdat;
}

sub slurp {
    # read the complete content of a file
    my ($file) = @_;
    local $/ = undef;
    open F, $file or die "Cannot read $file: $!\n";
    my $data = <F>;
    close F;
    return $data;
}
