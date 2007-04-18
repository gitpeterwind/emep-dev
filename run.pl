#!/usr/bin/perl

#Common script for Njord and Snykov. Choose $NJORD=1 or $SNYKOV=1

#___________________________________________________________________
#Snykov queue commands

#Queue system commands start with #PBS (these are not comments!)
# lnodes= number of nodes, ppn=processor per node (max4)
#PBS -lnodes=32
# wall time limit of run 
#PBS -lwalltime=01:00:00
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
# @ node = 2
# @ tasks_per_node = 16
# @ node_usage = not_shared
# @ network.mpi = sn_all,shared,us
# @ resources        = ConsumableCpus(1) ConsumableMemory(300 mb)
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

require "flush.pl";

#Choose one machine
$NJORD=0; #1 if njord is used
$SNYKOV=1; #1 if snykov (snowstorm) is used


$MAKE = "gmake --makefile=Makefile_snow";
  $MAKE = "make -f Makefile_njord" if $NJORD==1 ;
  die "Must choose SNYKOV **or** NJORD!\n" unless $NJORD+$SNYKOV==1;

$SR=0;     # Set to 1 if source-receptor calculation
 
# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables 
#      are explained below, and derived variables set later.-

$year = "2003";
( $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

$iyr_trend = $year;  
$iyr_trend = "2020" if $SR ;  # 2020 assumed for SR runs here

print "Year is $yy YEAR $year Trend year $iyr_trend\n";


#---  User-specific directories (changeable)

$PETER       = "mifapw/emep";      
$DAVE      = "mifads";      
$JOFFEN      = "mifajej";
$HILDE       = "mifahf";
$SVETLANA    = "mifast";
$HEIKO      = "mifahik";      
$ANNA      = "mifaab";      
$MICHAEL      = "mifaelg";      
$SEMEENA     = "mifasv";      
$TAREQ    = "mifatarh";      


$USER        =  $DAVE ;      
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

$DATA_LOCAL    = "$DataDir/EMEP";    # Grid specific data


my $HEMIS = 0;   #Set to 1 for Hemispheric run. Not possible yet
my $PM_ADDED     = 0;  # Adds in PM emissions from NOx inventory scaling
my $AFRICA_ADDED = 0;  # Adds in African emissions for y=1..11

$OZONE = "1"; 
$ACID = "0";     # Specify model type here, and check:
  die "Must choose ACID or OZONE" unless $OZONE+$ACID==1;


# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

if ( $OZONE ) {
    @emislist = qw ( sox nox nh3 co voc pm25 pmco ); 
    $testv       = "rv2_6_6";
    
} elsif ( $ACID ) {
    die "ACID not yet tested \n";	    
}


#User directories
$ProgDir     = "$HOMEROOT/$USER/Unify/Unimod.$testv";   # input of source-code
$WORKDIR     = "$WORKROOT/$USER/$testv.$year";    # working and result directory
$MyDataDir   = "$HOMEROOT/$USER/Unify/MyData";    # for each user's private input

#ds check: and change
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"  
#  unless $ENV{PWD} =~ /Unimod.$testv.$year/;


$Split       = "BASE_MAR2004" ;       
$NOxSplit       = "2000" ;               # Have CLE2020, MFR2020, 2000       
#$NOxSplit    = "CLE2020_ver2";    # NOx splits
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
$EMIS_INP = "$HOMEROOT/$HEIKO/Emissions/ModelInputModruns/Modrun06";
$emisdir = "$EMIS_INP/2006_emis2010_BL-E_V7";
$emisdir = "$EMIS_INP/2006-Trend${year}-V7";
$emisdir = "$EMIS_INP/2006-Trend2004-V7" if $year > 2004;
$pm_emisdir = $emisdir;
$pm_emisdir = "$EMIS_INP/2006-Trend2000-V7"  if $year < 2000;
 
$femis       = "$DataDir/femis.dat";      # emission control file


# Specify small domain if required. 
#                 x0   x1  y0   y1

@largedomain = (   1, 170,  1, 133 ) ;
##smalldomain = (  20, 167,  1, 122 ) ;    # OSPAR/HELCOM domain
##smalldomain = (  18, 169,  7, 124 ) ;     # OSPAR/HELCOM domain+border-south
@smalldomain = (  36, 167, 12, 122 ) ;    # EMEP domain
##@smalldomain = (  116, 167,  80, 122 ) ;      # (changeable)
##@smalldomain = @largedomain ;     # If you want to run for the whole domain, 
		# simply uncomment this  - REMEMBER NEED ##@, not just one #


$RESET        = 0   ;  # usually 0 (false) is ok, but set to 1 for full restart
$COMPILE_ONLY = 0   ;  # usually 0 (false) is ok, but set to 1 for compile-only
$INTERACTIVE  = 0   ;  # usually 0 (false), but set to 1 to make program stop
# just before execution - so code can be run interactivel.

$NDX   = 8;           # Processors in x-direction
$NDY   = 4;           # Processors in y-direction
if ( $INTERACTIVE ) { $NDX = $NDY = 1 };


@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

$mm1   =  "07";       # first month, use 2-digits!
$mm2   =  "07";       # last month, use 2-digits!
$NTERM_CALC =  calc_nterm($mm1,$mm2);

$NTERM =   $NTERM_CALC;    # sets NTERM for whole time-period
# -- or --
$NTERM = 4;       # for testing, simply reset here

print "NTERM_CALC = $NTERM_CALC, Used NTERM = $NTERM\n";

# <---------- end of normal use section ---------------------->
# <---------- start of SR   use section ---------------------->

#SR additions:
# Define all countries and nums here: Heiko's lib might be faster:
%country_nums = ( AT => 2, BE => 3, DK => 6, FI => 7, FR => 8,
                  DE => 10, GR => 11, IE => 14, IT => 15,
                  NL => 17 , PT => 20, ES => 22, SE => 23,UK => 27,
                  HU => 12, PL => 19, CY => 55, CZ => 46, EE => 43,
                  LT => 45, LV => 44, MA => 57, SK => 47, SI => 48 );
# EU countries:
@eu = qw ( AT BE DK FI FR DE GR IE IT NL PT ES SE UK );
@euaccess = qw ( HU PL CY CZ EE LT LV MA SK SI );
@eu25 = ( @eu, @euaccess );
@countries = qw ( IT ) ;  # List of wanted countries this run
# EU25 is "special"

$base        = "Base";

@polls       = qw ( BASE );  # NP, BASE 
    


# <---------- end of user-changeable section ----------------->
#               (normally, that is...)



#--- Verify data directories
foreach $d (  $WORKDIR, $DATA_LOCAL, $DataDir,  $ProgDir) {
    unless ( -d "$d" &&  -x _ && -r _ ) {
	die "*** ERROR *** directory $d not accessible. Exiting.\n";
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
    $smalldomain[2] < $largedomain[2] ||  $smalldomain[3] > $largedomain[3] );


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


if ( $RESET ) { ########## Recompile everything!
    
    
    # Set values for domain size in Par_ml.f90 : 
	open(EULPAR,"<Par_ml.pat") or die "No Par_ml.pat file!!\n";
	open(EULOUT,">Par_ml.f90");
	print "changing domain, nproc in Par_ml.pat to Par_ml.f90 \n";
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

die "Done. COMPILE ONLY\n" if  $COMPILE_ONLY;  ## exit after make ##

    
@list_of_files = ();   # Keep list of data-files



########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################
########################### START OF SCENARIO RUNS  ##########################


foreach $scenario ( @runs ) {
    print "STARTING RUN $scenario \n";
    
    $runlabel1    = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
    $runlabel2    = "${testv}_${scenario}_$year_$iyr_trend";   # NO SPACES! LONG (written into CDF files)
    
    my $RESDIR = "$WORKDIR/$scenario";
    system("mkdir -p $RESDIR");
    
    chdir "$RESDIR";   ############ ------ Change to RESDIR
    

for ($nnn = 1, $mm = $mm1; $mm <= $mm2; $mm++) {
    
    # Assign met files to fil001...etc.
    for ($n = 1; $n <= $month_days[$mm]; $n++) {
	$nnn = metlink($n, $nnn, $mm);
	last if $nnn > $NTERM;
    }
}

$mmlast = $mm2 + 1;
$yylast = $year;
if ( $mmlast > 12 && $NTERM > 200 ) { # Crude check that we aren't testing with NTERM=5
    $yylast = $yylast + 1;
    $mmlast = 1;
}
$old = sprintf "$MetDir/f00.%04d%02d01", $yylast, $mmlast;
$new = sprintf "fil%04d", $nnn;
mylink( "LAST RECORD SET: ", $old,$new ) ;


#=================== INPUT FILES =========================================
# ToDo Change noxsplit.default to defaults, as with voc (also in Unimod)
#   AFRICA, PM_ADDED --- fix elsewhere.

    %ifile   = ();   # List of input data-files

# First, emission files are labelled e.g. gridSOx, whiuch we assign to
# emislist.sox to ensure compatability with the names (sox,...) used
# in the model.

%gridmap = ( "co" => "CO", "nh3" => "NH3", "voc" => "NMVOC", "sox" => "SOx",
	 "nox" => "NOx" , "pm10" => "PM10", "pm25" => "PM25", "pmco" => "PMco" ) ;

    foreach $poll  ( @emislist  ) {
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
      $sondes="$DATA_LOCAL/sondes.dat";
      $sondes="$DATA_LOCAL/sondes.SR" if $SR;
    $ifile{"$sondes"} = "sondes.dat";
    $ifile{"$DATA_LOCAL/sites.dat"} = "sites.dat";
    $ifile{"$DATA_LOCAL/forests.mar2004b"} = "forest.dat";
    $ifile{"$DataDir/amilt42-nox.dat"} = "ancatmil.dat";#RENAME TO AIRCARAFT?!
    
# Seasonal stuff  ----    Can't we improve this? e.g. every month?
%seasons = ( "jan" => "01", "apr" => "02", "jul" => "03" , "oct" => "04") ;

    foreach $s ( keys(%seasons) ) {
	$ifile{"$DataDir/a${s}t42-nox.dat"} = "ancat$seasons{$s}.dat";
	$ifile{"$DataDir/jclear.$s"} = "jclear$seasons{$s}.dat";
	$ifile{"$DataDir/jcl1.$s"} = "jcl1km$seasons{$s}.dat";
	$ifile{"$DataDir/jcl3.$s"} = "jcl3km$seasons{$s}.dat";
    } 
    
    $ifile{"$DATA_LOCAL/rough.170"} = "rough.170"; # Roughness length;
    $ifile{"$DATA_LOCAL/landuse.JUN06"} = "landuse.JUN06";
    $ifile{"$DATA_LOCAL/Volcanoes.dat"} = "Volcanoes.dat";

    $ifile{"$DataDir/JUN06_gfac1.dat"} = "JUN06_gfac1.dat";
    $ifile{"$DataDir/JUN06_gfac2.dat"} = "JUN06_gfac2.dat";
    $ifile{"$DataDir/JUN06_biomass.dat"} = "JUN06_biomass.dat";

    foreach my $old ( sort keys %ifile ) {  # CHECK and LINK
	if ( -r $old ) {
		$new =  $ifile{$old};
       		mylink( "Inputs: ", $old,$new ) ;
	} else {
        	print "Missing Input $old !!!\n";
		die "ERROR: Missing OLD $old\n" unless $old =~ /special/;
	}
    }
#=================== INPUT FILES =========================================

    
# FIX later - was the only emission control thingy....
    @exclus  = (9 ); #  NBOUND
    
    
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
    
    &flush(STDOUT);
    
# Link execultable also, since gridur is funny about these
# things....
    
    $LPROG = "prog.exe";
    #mylink( "PROGRAM!!  ", $PROGRAM,$LPROG) ;
    system( "cp $PROGRAM $LPROG") ;
    
# Write out list of linked files to a shell-script, useful in case the program
# hangs or crashes:
    
    open(RMF,">Remove.sh");
    foreach $f ( @list_of_files ) { print RMF "rm $f \n" };
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
	

	my $PRERUN = "";
	$PRERUN = "scampiexec " if $SNYKOV;
        open (PROG, "| $PRERUN ./$LPROG") || 
		die "Unable to execute $LPROG. Exiting.\\n" ;
        close(PROG);
	    
	
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
        print "\n  Unkown Problem!! \n";
    }

#move RunLog 
    system("mv RunLog.out  ${runlabel1}_RunLog");
    system("echo Emission units: Gg/year        >> ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    system("echo Emissions: $emisdir            >> ${runlabel1}_RunLog");
    system("echo Version: $testv                >> ${runlabel1}_RunLog");
    system("echo Domain x0 $dom_x0 y0 $dom_y0 wx $dom_wx wy $dom_wy  >> ${runlabel1}_RunLog");
    system("echo Processors $NDX $NDY           >> ${runlabel1}_RunLog");
    system("echo Added? PM $PM_ADDED  Africa $AFRICA_ADDED >> ${runlabel1}_RunLog");
    system("echo SR?  $SR                       >> ${runlabel1}_RunLog");
    system("echo iyr_trend: $iyr_trend                 >> ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    system("echo femis: femis.$scenario         >> ${runlabel1}_RunLog");
    system("cat femis.dat >> ${runlabel1}_RunLog");
    system("echo ------------------------------ >> ${runlabel1}_RunLog");
    
    
#clean up work directories and links   

    unlink ( @list_of_files );

#tar sites and sondes. Use sondes to check as these are produced les frequently.
    my $last_sondes = sprintf  "sondes.%02d%02d", $mm2, $yy;
    print "LOOKING FOR LAST SITES $last_sondes\n";
    if ( -r $last_sondes ) {
	print "FOUND LAST sondes $last_sondes\n";
	system("tar cvf ${runlabel1}.sites sites.*");
	system(	"tar cvf ${runlabel1}.sondes sondes.*");
    }
    
    
################################## END OF SCENARIO RUNS ######################
}  ############################### END OF SCENARIO RUNS ######################
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

sub metlink {  #---- meteorological data
    local ($dd,$nnn,$mm) = ($_[0], $_[1], $_[2]);

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

