#!/usr/bin/perl

#Common script for Njord, Stallo.
#Choose $VILJE=1 or $STALLO=1
#___________________________________________________________________

#Queue system commands start with #SBATCH for Stallo and #PBS for Vilje (these are not comments!)
# to activate take out one # from ##
#___________________________________________________________________
#Stallo SLURM queue commands
#submit with:  sbatch run.pl
#for express queue (max 4 hours):  sbatch --qos=devel run.pl
#Queue system commands start with #SBATCH 
#SBATCH -A nn2890k
#you can change number of nodes. Total CPU = nodes*ntasks-per-node
#SBATCH --nodes=2 
#SBATCH --ntasks-per-node=20
##SBATCH --ntasks=40 Do not use for now! (bug on Stallo)
#SBATCH --mem=32000
#SBATCH --time=4:0:0
#activate the following line for runs which last longer than 48 hours AND use more than one node
##SBATCH --partition=multinode 
#SBATCH --job-name=emep
#SBATCH --output=run.%j.out

# Vilje:
#   select= number of nodes, ncpus=number of threads per node to reserve, 
# to activate take out one # from ##
#   mpiprocs=number of MPI threads per node. select=number of nodes
# Q. 32 or 16 threads? There are only 16 cores per node, so if asking for 32 threads,
# the 2 threads per core run on the same core and try to utilize different part
# of the core (like one thread fetch from memory, while the other thread does
# a multiplication). This is not very efficient, therefore it will not run
#  much faster than if you use select=4:ncpus=16:mpiprocs=16) 
#
#PBS -l select=4:ncpus=16:mpiprocs=16 -v MPI_MSGS_MAX=2097152,MPI_BUFS_PER_PROC=2048
#
# Wall time limit of run
#PBS -lwalltime=07:20:00
# Make results readable for others:
#PBS -W umask=0022
# Account for billing
#PBS -A nn2890k
##PBS -A mifa01kl
# Multiple tasks for paralel SR runs (one task per country)
##PBS -t 1-56
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
######################################################################
#Tips:            SNYKOV or STALLO  NJORD            ABEL?
#  submit job     qsub run.pl       llsubmit run.pl  sbatch run.sh
#  queue status   qstat -u $USER    llq              squeue -u $USER
#  job status                       checkjob 3456    checkjob 3456
#  kill job       qdel 3456         llcancel 3456    scancel 3456
#  submit multi-task SR job (one task per country)
#                 qsub -t 1-56 run.pl                arrayrun 1-56 run.sh
######################################################################

use 5.6.0;
use feature qw{ switch };
use strict;
use warnings;
use File::Copy qw();
use File::Compare;

$| = 1; # autoflush STDOUT

#Choose one machine
#my $VILJE=0;  #1 if Ve or Vilje is used
#my $STALLO=1; #1 if stallo is used
my ( $STALLO, $VILJE ) = (0) x 2;
foreach my $key  (qw ( PBS_O_HOST HOSTNAME PBS_SERVER MACHINE SLURM_SUBMIT_HOST)) {
 next unless defined $ENV{$key};
 $STALLO = 1 if $ENV{$key} =~/stallo/;
 $VILJE  = 1 if $ENV{$key} =~/vilje/;
 #my $val = $ENV{$key};
 #print "KEY $key Val $val  STALLO?? $STALLO VILJE $VILJE\n";
}
print "HOST DETECT: V$VILJE S$STALLO \n";

# -j4 parallel make with 4 threads
my @MAKE = ("gmake", "-j8", "MACHINE=snow");
   @MAKE = ( "make", "-j8", "MACHINE=vilje")  if $VILJE==1 ;
   @MAKE = ( "make", "-j8", "MACHINE=stallo") if $STALLO==1 ;
die "Must choose STALLO **or** VILJE !\n"
  unless $STALLO+$VILJE==1;

# Setup of code, domain,  outputs, chemistry
# $Chem: EmChem09soa, EmChem09, CRI_v2_R5
# $GRIDs:EECCA, EMEP, TNO7, TNO14, TNO28, TNO56, MACC02, MACC14 or GLOBAL
# MAKEMODE - be careful!
# By default $MAKEMODE is set to EMEP, which implies EmChem09soa.
# $MAKEMODE: EMEP,EMEP2011,eEMEP,MACC,MACC-EVA2011,... see Makefile
# $MAKEMODE=0; # Avoid model re-compilation (skips call to GenChem)
#
# You could set $MAKEMODE to "EmChem09". For instance:
# my ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) =
#    ("2674","EmChem09","EMEPSTD","EMEPSTD","EECCA","EmChem09");
#
my ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("rv4_6gamma"   ,"EmChem09soa","EMEPSTD","EMEPSTD","EECCA","EMEP");
#  ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("test"    ,"EmChem09"   ,"EMEPSTD","EMEPSTD","EECCA",0);
#  ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("testcri2","CRI_v2_R5"  ,"CRITEST","EMEPSTD","EECCA",0);
#eg ($testv,$Chem,$exp_name,$GRID,$MAKEMODE) = ("tests","EmChem09","TESTS","RCA","EmChem09");
 ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("rv4_11_3","EmChem09soa","EMEPSTD","EMEPSTD","EECCA",0);
#($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("3074","EmChem09soa","EMEPGLOB","EMEPSTD","GLOBAL",0);
 ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("emep-dev","EmChem09soa","EMEPSTD","EMEPSTD","EECCA",0);
 ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("rv4_17box","EmChem16x","A2018","EMEPSTD","EECCA",0);
 ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("rv4_17aug13","EmChem16x","A2018e","EMEPSTD","EECCA",0);

my %BENCHMARK;
# OpenSource 2008
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2005,emis=>"Modrun07/OpenSourceEmis"           ,chem=>"EmChem03");
# Dave's preference for EMEP:
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2006,emis=>"Modrun10/EMEP_trend_2000-2008/2006",chem=>"EmChem09");
# EECCA Default:
# %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun11/EMEP_trend_2000-2009/2008",chem=>"EmChem09soa",make=>"EMEP");
# Status Runs:
#  %BENCHMARK = (grid=>"EECCA" ,year=>2007,emis=>"Modrun09/2009-Trend2007-CEIP") ;
#  %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun10/2010-Trend2008_CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2009,emis=>"Modrun11/2011-Trend2009-CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2010,emis=>"Modrun12/2012-Trend2010-CEIP",chem=>"EmChem09soa",make=>"EMEP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2011,emis=>"Modrun13/2013-Trend2011-CEIP",chem=>"EmChem09soa",make=>"EMEP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2012,emis=>"Modrun14/2014-Trend2012-CEIP",chem=>"EmChem09soa",make=>"EMEP");
# Alternative domains:
#  %BENCHMARK = (grid=>"TNO28" ,year=>2008,emis=>"emis_TNO28"         );
#  %BENCHMARK = (grid=>"MACC02",year=>2008,emis=>"2008_emis_EMEP_MACC") ;
if (%BENCHMARK) {
  $BENCHMARK{'archive'} = 1;                   # save summary info in $DataDir
  $BENCHMARK{'debug'} = $BENCHMARK{'archive'}; # chech if all debug flags are .false.
# Default setting, if not previously specified
  $BENCHMARK{'chem'}  = "EmChem09soa"
    unless $BENCHMARK{'chem'};  # chemical mecanism, e.g. OpenSource 2008
  $BENCHMARK{'make'}  = ($BENCHMARK{'chem'} eq "EmChem09soa")?"EMEP":"all"
    unless $BENCHMARK{'make'};  # make target, e.g. Status 2010
  $Chem     = $BENCHMARK{'chem'};
  $MAKEMODE = $BENCHMARK{'make'};
  $GRID     = $BENCHMARK{'grid'};
}

my $INERIS_FACS=0;  # Used for timefactors, and e,g TNOxx tests
my $SR= 0;     # Set to 1 if source-receptor calculation
               # check also variables in package EMEP::Sr below!!
 $MAKEMODE="SR-$MAKEMODE" if($MAKEMODE and $SR);

# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables
#      are explained below, and derived variables set later.-

my $year = "2012";
   $year = $BENCHMARK{"year"} if %BENCHMARK;
( my $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

my $iyr_trend = $year;
#$iyr_trend = "1990" ;  #  RCA ECLAIRE

print "Year is $yy YEAR $year Trend year $iyr_trend\n";


#---  User-specific directories (changeable)

my $PETER      = "mifapw";
my $DAVE       = "mifads";
my $JOFFEN     = "mifajej";
my $HILDE      = "mifahf";
my $SVETLANA   = "mifast";
my $HEIKO      = "mifahik";
my $ANNA       = "mifaab";
my $MICHAEL    = "michaelg";
my $SEMEENA    = "mifasv";
my $AGNES      = "nyiri";
my $ALVARO     = "alvarov";
my $ROBERT     = "mifarb";
my $HALDIS     = "mifahb";
my $BIRTHE     = "birthems";
my $FORCAST    = "forecast";

my $USER = $ENV{"USER"};
print "USER = $USER\n";


# hb NH3Emis
my $NH3EMIS_VAR = 0; # set to 1 if new temp NH3.

my $MONTHLY_EMIS = ($GRID eq "GLOBAL"); #Switch off if only annual used
   $MONTHLY_EMIS = 0 if ($exp_name =~ /ECLAIRE/);

my ($HOMEROOT, $WORKROOT, $MetDir);
our $DataDir;
if ($STALLO) {
  $HOMEROOT = "/home";
  $WORKROOT = "/global/work";
  $DataDir  = "$WORKROOT/$PETER/emep/Data";
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year";
  $MetDir   = "$DataDir/$GRID/metdata_H20/$year" unless -d $MetDir;
  $MetDir   = "$DataDir/$GRID/metdata/$year"     unless -d $MetDir;
  $MetDir   = "$WORKROOT/$PETER/emep/ClimData/$year"  if ($GRID eq "RCA" );
} else { #Ve or Vilje
  $HOMEROOT = "/home/metno";
  $WORKROOT = "$HOMEROOT/$USER/work";
  $DataDir  = "$HOMEROOT/mifapw/work/Data";
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year";
}
#$MetDir   = "/hard/code/path/to/$GRID/metdata/$year/if/necessary";
die "Missing MetDir='$MetDir'\n" unless -d date2str($year."0101",$MetDir);

# DataDir    = Main general Data directory
my $DATA_LOCAL = "$DataDir/$GRID";   # Grid specific data , EMEP, EECCA, GLOBAL

# Project-specific data:
my $ProjDataDir = "";          # Change for specific project data
#  $ProjDataDir = "$WORKROOT/mifads/Data/inputs_projects/eclaire_Jun2013"; #e.g.

# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

my $VBS   = 0;

#User directories
my $ProgDir  = "$HOMEROOT/$USER/emep-mscw";   # input of source-code
   $ProgDir  = "$HOMEROOT/$USER/emep-mscw/$testv";   # input of source-code for testv
my $ChemDir  = "$ProgDir/ZCM_$Chem";
my $Specials = "specials";  # default
#$Specials = "TSAP_Jul2012";  # used for TSAP runs in July 2012
# Check:
die "No ProgDir! $ProgDir\n" unless -d $ProgDir;
die "No ChemDir! $ChemDir\n" unless -d $ChemDir;


#---- emislist --------------------------------------------------------
#A2018 skip: open(EMIS,"<$ProgDir/CM_emislist.csv") or die "Need CM_emislist.cvs file!\n";
#  my @emislist = split(/,/,<EMIS>);
#  print "EMISLIST ", join(" ", @emislist ), "\n";
#close(EMIS);
#A2018 skip
#----  chem packages  (e.g. EmChembase PMmass ) -----------------------
#open(CHEM,"<$ProgDir/CM_chempackages.txt") or die "Need CM_emislist.cvs file!\n";
#  my @packages = <CHEM> or die "Need CM_chempackage.txt!\n" ;
#  print "CHEM packages:\n @packages\n";
#close(CHEM);
#----------------------------------------------------------------------

# Check that the code directory has the chem files we want:
# Now use mk.GenChem, CM_ files erased from ZCM_ directories
#die "Mis-Match chemistry, emepctm.$testv Chem: $Chem" if
#  ( File::Compare::compare( "$ProgDir/CM_ChemSpecs_ml.f90" , "$ChemDir/CM_ChemSpecs_ml.f90"));

my $WORKDIR = "$WORKROOT/$USER/$testv.$year";  # working and result directory
   $WORKDIR = "$WORKROOT/$testv.$year" if($WORKROOT =~ /$USER/);
   $WORKDIR =~ s|$testv.$year|Benchmark/$GRID.$year|g if (%BENCHMARK);
my $MyDataDir="$HOMEROOT/$USER/Unify/MyData"; # for each user's private input
my $SoilDir = "$DATA_LOCAL/dust_input";       # Saharan BIC
   $SoilDir = 0 unless -d "$SoilDir/BC_DUST/2000";

# TEST! Road dust NOTE! The road dust coddumpe may not be working properly yet! Not tested enough!
my $RoadDir = "$HOMEROOT/$ROBERT/Unify/MyData/TNO_traffic" ;
   $RoadDir = "$HOMEROOT/$AGNES/MyData/TNO_traffic" if $VILJE;
   $RoadDir = 0 unless -d $RoadDir;

#ds check: and change
chdir "$ProgDir";
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"
#  unless $ENV{PWD} =~ /emepctm.$testv.$year/;
print "TESTING ENV:", $ENV{PWD}, "\n";

my $version     = "emepctm" ;
my $PROGRAM     = "$ProgDir/$version";        # programme
my $subv        = "$testv" ;                  # sub-version (to track changes)

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later.
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR
my @runs     = ( $scenario );


my $RESET        = 0; # usually 0 (false) is ok, but set to 1 for full restart
my $COMPILE_ONLY = 0; # usually 0 (false) is ok, but set to 1 for compile-only
my $DRY_RUN      = 0; # Test script without running model (but compiling)
my $KEEP_LINKS   = 0; # Keep @list_of_files after run
#$KEEP_LINKS=not $BENCHMARK{'archive'} if(%BENCHMARK);

if(%BENCHMARK and $BENCHMARK{'debug'}){
  die "No debug flags for benchmarks!"
  if(system("grep -Hnie 'logical.*DEBUG.*=\ *.TRUE.' $ProgDir/*.f90")==0) or
    (system("grep -Hnie 'DEBUG.*=\ *.TRUE.' $ProgDir/ModelConstants_ml.f90")==0);
}

if($ENV{PBS_NODEFILE}){
  $_ =  `wc -l $ENV{PBS_NODEFILE}`;
  my $RUN_NPROC;
  ($RUN_NPROC,undef) = split;
  print "Qsub required: $RUN_NPROC processors\n";
}else{
  print "skip nodefile check on interactive runs\n";
}

my @month_days=(0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

#Only 360 days in HIRHAM metdata. We ignore leaps
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,24) if $GRID eq "HIRHAM";

my $mm1 ="06";      # first month, use 2-digits!
my $mm2 ="06";      # last month, use 2-digits!
my $dd1 =  1;       # Start day, usually 1
my $dd2 =  2;       # End day (can be too large; will be limited to max number of days in the month)
                    # put dd2=0 for 1 timestep run/test.
# Allways runn full year on benchmark mode
($mm1,$mm2,$dd1,$dd2)=("01","12",1,31) if (%BENCHMARK);
$dd1=($dd1>$month_days[$mm1])?$month_days[$mm1]:$dd1;
$dd2=($dd2>$month_days[$mm2])?$month_days[$mm2]:$dd2;

# <---------- end of normal use section ---------------------->
# <---------- end of user-changeable section ----------------->
#               (normally, that is...)

if ($SR) {
    print "SR is true\n";
    @runs = EMEP::Sr::initRuns();
}

if (%BENCHMARK){
  print "BENCHMARK is true: $GRID $year $testv $Chem\n";
  @runs = ( "$testv-$Chem" );
}


#--- Verify data directories
mkdir_p($WORKDIR);
foreach my $d (  $WORKDIR, $DataDir,  $ProgDir) {
    unless ( -d "$d" &&  -x _ && -r _ ) {
        die "*** ERROR *** directory $d not accessible. Exiting.\n";
    }
}


#--- Other eulmod configs

#
# Check that we have an existing prog dir:
die "Wrong ProgDir: $ProgDir \n" unless -d $ProgDir;



#--- Calendar stuff
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);


#--- adjust for leap year
$month_days[2] += leap_year($year);

# --- Start new compilations if needed

# --- We check if the already-compiled version of $PROGRAM is the same
#     as the one we are asking for. The old configuration should have
#     been stored in Make.log.  Read in this and split terms into
#     : model, emis


chdir "$ProgDir";

#-- generate Makefile each time, to avoid forgetting changed "pat" file!

if($RESET) { ########## Recompile everything!
  # For now, we simply recompile everything!
  system(@MAKE, "clean");
  if ($MAKEMODE) {
    system(@MAKE, "clean", "$MAKEMODE") if ($MAKEMODE=~/3DVar/);
    system(@MAKE, "$MAKEMODE");
  } else {
    system(@MAKE, "depend");
    system(@MAKE, "all");
  }
}
system "pwd";
print "Check last files modified:\n";
system "ls -lht --time-style=long-iso -I\*{~,.o,.mod} | head -6 ";

#to be sure that we don't use an old version (recommended while developing)
#unlink($PROGRAM);

if ($SR or ($USER eq $FORCAST) or
    defined($ENV{'PBS_ARRAY_INDEX'}) or 
    defined($ENV{'PBS_ARRAYID'}) or 
    defined($ENV{'TASK_ID'})){
  #No recompile SR or ARRAY jobs
} elsif ($MAKEMODE) {
  system(@MAKE, "$MAKEMODE") == 0 or die "@MAKE $MAKEMODE failed";
} else {
  system (@MAKE, "depend") ;
  system (@MAKE, "all") == 0 or die "@MAKE all failed";
}

die "Done. COMPILE ONLY\n" if  $COMPILE_ONLY;  ## exit after make ##


my @list_of_files = ();   # Keep list of data-files

########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################

foreach my $scenflag ( @runs ) {
  if ($SR) {
    $scenario = EMEP::Sr::getScenario($scenflag);
  } elsif (%BENCHMARK){
    $scenario = "BM_$scenflag";
  } else {
    $scenario = $scenflag;
  }
  print "STARTING RUN $scenario \n";

  my $runlabel1 = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
  my $svn = "$ProgDir/.version";
  chop($svn = (-e $svn)?qx(cat -s $svn):qx(svnversion -n));
  my $runlabel2 = "$testv\_$Chem\_svn$svn\_$scenario\_$year\_Trend$iyr_trend";   # NO SPACES! LONG (written into CDF files)

  my $RESDIR = "$WORKDIR/$scenario";
     $RESDIR = "$WORKDIR/$scenario.$iyr_trend" if ($GRID eq "RCA");
  mkdir_p($RESDIR);

  chdir $RESDIR;   ############ ------ Change to RESDIR
  print "Working in directory: $RESDIR\n";

#To use FastJ some data files are required. Could be moved elsewhere 
 my $FASTJ = 0;
 if($FASTJ) { 
  my $old = "$DataDir/fastj/FJX_j2j.dat";
  my $new = "FJX_j2j.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_scat-aer.dat";
  $new = "FJX_scat-aer.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_scat-cld.dat";
  $new = "FJX_scat-cld.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_scat-UMa.dat";
  $new = "FJX_scat-UMa.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_spec-200nm-2013c.dat";
  $new = "FJX_spec-200nm-2013c.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_spec.dat";
  $new = "FJX_spec.dat";
  mylink( "Linking:", $old, $new);
  $old = "$DataDir/fastj/FJX_spec-std-2013c.dat";
  $new = "FJX_spec-std-2013c.dat";
  mylink( "Linking:", $old, $new);
 }
#=================== INPUT FILES =========================================
# ToDo Change noxsplit.default to defaults, as with voc (also in emepctm)


  EMEP::Sr::generate_updated_femis(@$scenflag) if ($SR);
#=================== INPUT FILES =========================================
# FIX later - was the only emission control thingy....
  my @exclus  = (9 ); #  NBOUND
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
#------------      Run model ------------------------------------------
  print "\n";

# compile/copy executable to $RESDIR/
  my $LPROG = "emepctm";
  if($MAKEMODE){  # compile/link program into $RESDIR/
    system(@MAKE,"$MAKEMODE","-C$ProgDir/",
           "ARCHIVE=yes","BINDIR=$RESDIR/") == 0
      or die "@MAKE $MAKEMODE -C$ProgDir/ failed"; 
  }else{          # copy program to $RESDIR/
    cp ($PROGRAM, $LPROG) or die "cannot copy $PROGRAM to $LPROG\n";
  }
  push(@list_of_files , $LPROG);    # For later deletion

# Write out list of linked files to a shell-script, useful in case the program
# hangs or crashes:

  open(RMF,">Remove.sh");
  foreach my $f ( @list_of_files ) { print RMF "rm $f \n"; }
  print RMF "rm $LPROG\n";   # Delete executable also
  close(RMF);

  my ($startdate,$enddate)=("$year-$mm1-$dd1","$year-$mm2-$dd2");
     $enddate=date2str($startdate." 1 day ago","%F") unless $dd2;
  $startdate=date2str("$startdate","%Y%m%d");
  $enddate  =date2str("$enddate"  ,"%Y%m%d");
# check if met-files exist
#  for (my $d="$startdate";$d<=$enddate;$d=date2str($d." 1 day","%Y%m%d")) {
#    my $f=date2str($d,$METformat);
#    die "METFILE not found:\n\t$f\n" unless -e $f;
#  }

# namelist with input parameters (to be read by emep_Main.f90 and other modules)
  my $nml="";
  my %h=(%BENCHMARK,MetDir=>"$MetDir");
  if(%BENCHMARK){
    # read nml template file
    $nml=EMEP::Sr::slurp("$ProgDir/config_BM-$GRID.nml");
    # fill in variables on the template file with corresponding $hash{key}
    %h=(%h,'runlabel1'=>"$runlabel1",'runlabel2'=>"$runlabel2");
  } else {
    $nml="&INPUT_PARA\n"
        ."  GRID      = '$GRID',\n"
        ."  iyr_trend = $iyr_trend,\n"
        ."  runlabel1 = '$runlabel1',\n"
        ."  runlabel2 = '$runlabel2',\n"
        ."  startdate = ".date2str($startdate ,"%Y,%m,%d,000000,\n")
        ."  enddate   = ".date2str($enddate   ,"%Y,%m,%d,000024,\n")
#       ."  meteo     = '$METformat',\n" #moved to config
        ."&end\n";
    # NML namelist options.
    foreach my $f ("config_$exp_name.nml","config_Outputs_$outputs.nml") {
      $nml.=EMEP::Sr::slurp("$ProgDir/$f");
    }
    # fill in variables on the template file with corresponding $hash{key}
  }
  # fill in variables on the template file with corresponding $hash{key}
  foreach my $k (keys %h) {
    $nml=~s:\$$k:$h{$k}:g;  # replace keyword ('$key') by its value $h{'key'}
    my $dk=(-e $h{$k})?`date -r $h{$k} +"\#%F %R %Z"`:"\n";
    $nml=~s:\#$k:$dk:g;     # and ('#file') by its time stamp
  }
  # list mode setup variables
  $nml.="#". "-"x22 ." Model set-up ". "-"x22 ."\n";
  %h=('testv'=>"$testv",'Chem'=>"$Chem",'exp_name'=>"$exp_name",
      'outputs'=>"$outputs",'GRID'=>"$GRID",'MAKEMODE'=>"$MAKEMODE");
  foreach my $k (sort keys %h) {
    $nml.=sprintf "#  %-22s = '%s',\n",$k,$h{$k};
  }
  # list of linked files in nml format for compatibility with future nml-only versions
  $nml.="#". "-"x22 ." Linked files ". "-"x22 ."\n";

  $nml =~ s/(\s*\!.*|\s+$)//g;  # remove comments (!) and tailing spaces
  $nml =~ s/\#/\!/g;            # annotations (#) as comments (!)
  $nml =~ s/\s*\n+/\n/g;        # remove empty lines
  $nml.="\n";                   # restore newline at the end of last namelist
  open(TMP,">config_emep.nml");
  print TMP "$nml";
  close(TMP);

  foreach my $exclu ( @exclus) {
    print "starting $PROGRAM with\n".
    "EXCLU $exclu\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2\n".
    "startdate $startdate\nenddate $enddate\n";

    my $MPIRUN = "mpiexec";
    if ($STALLO) {
      $MPIRUN = "mpirun";
    } elsif ($USER eq $FORCAST) {
      $MPIRUN = "mpiexec_mpt"; # forecast user use special mpiexec_mpt!
    }
    if ($DRY_RUN) {
      print "DRY_RUN: not running '| $MPIRUN ./$LPROG'\n";
      system("ls -lht --time-style=long-iso *");
    } else {
      open (PROG, "| $MPIRUN ./$LPROG") ||
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
------------------------------
Emission units: Gg/year
------------------------------
Version: $testv
Chemical scheme: $Chem
SR?  $SR
iyr_trend: $iyr_trend
------------------------------
femis: femis.$scenario
------------------------------
EOT
  close RUNLOG;
  if ( -s "$runlabel1\_fullrun.nc") {
    system("cdo infov $runlabel1\_fullrun.nc > $runlabel1\_fullrun.infov");
  }

# BENCHMARK summary info
  if (%BENCHMARK and -s "$runlabel1\_fullrun.nc") {
    system("ls -lh --time-style=long-iso * > $runlabel1\_filelist");
    if ($BENCHMARK{'archive'}) {
      my @fileBM=( "$runlabel1\_fullrun.infov", "$runlabel1\_RunLog",
                   "$runlabel1\_filelist", "eulmod.res", "Timing.out" );
      my $dirBM="$DataDir/Benchmark/$GRID.$year/BM_$testv-$Chem.$USER";
      mkdir_p($dirBM);
      foreach my $f ( @fileBM ) { cp ("$f","$dirBM/$f") if -e "$f"; }
    }
  }

#clean up work directories and links
  if ($DRY_RUN or $SR){ # keep femis.dat
    @list_of_files = grep {$_ ne 'femis.dat'} @list_of_files;
  }
 unlink ( @list_of_files ) unless($KEEP_LINKS);


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

sub date2str { 
  my ($date,$str) = ("$_[0]","$_[1]");
  $date =~ s|,000000| 00:00:00|g;
  $date =~ s|,0000| 00:00|g;
  $date =~ s|,|-|g;
  $str =~ s|YYYY|%Y|g;
  $str =~ s|MM|%m|g;
  $str =~ s|DD|%d|g;
  chop($str = `date -d '$date' +'$str'`);
  return $str;
}

sub mylink {
  # links files from the original location (old) to
  # the new location (new) - generally the working directory.
  # Keeps track of all such linked files in list_of_files.
  my ($text,$old,$new) = ($_[0],$_[1],$_[2]);
  # remove new file if already exists and is symlink
  unlink $new if -l $new;
  symlink ($old,$new) || die "symlink $old $new failed : $!";
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

##############################################################
##### Stuff for Source receptor matrisses             ########
##############################################################
package EMEP::Sr;

my (%country_nums, %city_lonlat, @eu15, @euNew04, @eu25, @euNew06, @eu27, @sea, @noneu, @emep, @eecca, @eccomb);
our ($base, $Split, $NOxSplit, $rednflag, $redn, @countries, @polls);

INIT {
########################################
# Define all countries and nums here: ##
########################################
%country_nums = (
  AL =>  1,AT =>  2,BE =>  3,BG =>  4,FCS=>  5,DK =>  6,FI =>  7,FR =>  8,
  GR => 11,HU => 12,IS => 13,IE => 14,IT => 15,LU => 16,NL => 17,NO => 18,PL => 19,PT => 20,
  RO => 21,ES => 22,SE => 23,CH => 24,TR => 25,FSU=> 26,GB => 27,         REM=> 29,BAS=> 30,
  NOS=> 31,ATL=> 32,MED=> 33,BLS=> 34,NAT=> 35,RUO=> 36,RUP=> 37,RUA=> 38,BY => 39,UA => 40,
  MD => 41,RUR=> 42,EE => 43,LV => 44,LT => 45,CZ => 46,SK => 47,SI => 48,HR => 49,BA => 50,
  CS => 51,MK => 52,KZ => 53,GE => 54,CY => 55,AM => 56,MT => 57,ASI=> 58,LI => 59,DE => 60,
  RU => 61,MC => 62,NOA=> 63,EU => 64,US => 65,CA => 66,BIC=> 67,KG => 68,AZ => 69,ATX=> 70,
  RUX=> 71,RS => 72,ME => 73,RFE=> 74,KZE=> 75,UZ => 76,TM => 77,UZE=> 78,TME=> 79,CAS=> 80,
  TJ => 81,ARL=> 82,ARE=> 83,ASM=> 84,ASE=> 85,AOE=> 86,RFX=> 87,ASX=> 88,PAX=> 89,AOX=> 90,
  NAX=> 91,KZT=> 92,RUE=> 93,UZT=> 94,TMT=> 95,AST=> 96,
  BA2=>302,BA3=>303,BA4=>304,BA5=>305,BA6=>306,BA7=>307,BA8=>308,BA9=>309, # Baltic sep.  
  NS2=>312,NS3=>313,NS4=>314,NS5=>315,NS6=>316,NS7=>317,NS8=>318,NS9=>319, # N. Sea sep.
  AT2=>322,AT3=>323,AT4=>324,AT5=>325,AT6=>326,AT7=>327,AT8=>328,AT9=>329, # Atlant.sep.
  ME2=>332,ME3=>333,ME4=>334,ME5=>335,ME6=>336,ME7=>337,ME8=>338,ME9=>339, # Medit. sep.
  BL2=>342,BL3=>343,BL4=>344,BL5=>345,BL6=>346,BL7=>347,BL8=>348,BL9=>349, # Bl.Sea sep.
  ALL=>  0
);

# EU countries:
@eu15 = qw ( AT BE DK FI FR DE GR IE IT NL PT ES SE GB LU );
@euNew04 = qw ( HU PL CY CZ EE LT LV MT SK SI );
@eu25 = ( @eu15, @euNew04 );
@euNew06 = qw(BG RO);
@eu27 = (@eu25, @euNew06);
@sea = qw ( NOS ATL MED BAS BLS );
@noneu = qw ( NO CH IS );
@emep = qw ( RS ME BY BA HR TR UA MD MK GE AM AL AZ NOA ASI) ;
@eecca = qw ( KG KZ RU TJ );
@eccomb = qw ( RUE KZT UZT TMT AST );
########################################
# End of country definitions          ##
########################################


################################
#### start of SR parameters ####
################################
$base        = "CLE";
#$Split       = "CLE_MAR2004";     # Defualt (IER-based) VOC splits
#$NOxSplit    = "CLE2020_ver2";    # Default scenario (IER-based) VOC splits
$rednflag    = "P15";  # 15% reduction label
$redn        = "0.85"; # 15% reduction

# modify those to fill up your queues for SR effectively!!!
@countries  = (@eu27, @sea, @noneu, @emep, @eecca, @eccomb);
@polls       = qw (BASE N P A V S );  #  (any, all, at least 1)

# multiple tasks for paralel SR runs: one task per country
# Queue system:    POSIX             PBS           SLURM
foreach my $task ('PBS_ARRAY_INDEX','PBS_ARRAYID','TASK_ID') {
  @countries=($countries[$ENV{$task}-1]) if $ENV{$task};
}

if(defined($ENV{'CitySR'})){
# read City_Domain.dat file into %city
  my %city = ();
  my $fcity=slurp($ENV{'CitySR'});
  foreach my $line (split(/\n/,$fcity)){ # Split lines
    my ($cc,$ln1,$ln2,$lt1,$lt2)=split(" ",$line);
     $city{"$cc"}="$ln1 $ln2 $lt1 $lt2";
  }
# city pairs for femis.dat
  %city_lonlat = (
    ParBer=>["lonlat $city{Paris}",
             "lonlat $city{Berlin}"],
    LonPoV=>["lonlat $city{London}",
             "lonlat $city{PoValley}"],
    OslRuh=>["lonlat $city{Oslo}",
             "lonlat $city{Ruhrgebiet}"]
  );
  @countries=("ALL",(keys %city_lonlat));
  @polls    =("BASE","ALL");
# $rednflag ="CitySR";

# multiple tasks for paralel SR runs: one task per ciry pair
# Queue system:    POSIX             PBS           SLURM
  foreach my $key ('PBS_ARRAY_INDEX','PBS_ARRAYID','TASK_ID') {
    next unless defined $ENV{$key};
    if($ENV{$key} eq 0){
      @countries=($countries[0]);
      @polls    =("BASE");
    }else{
      @countries=($countries[$ENV{$key}-1]);
      @polls    =("ALL");
    }
  }
}
################################
#### end of SR parameters   ####
################################
}

sub initRuns {
  my @runs;
  foreach my $cc (@countries) {
    foreach my $poll (@polls) {
      push @runs, [$cc, $poll, $redn];
    }
    # run BASE only once (for exactly one cc)!!!
    @polls = grep {'BASE' ne $_} @polls;
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
  my ($cc,$pollut,$redn) = @_;
  # Initialise to 1.0:
  my($sox,$nox,$voc,$nh3,$testp,$co,$pm25,$pmco) = ("1.0")x8 ;
  given($pollut){
    when("AV"  ){ $voc = $nh3 = $redn; }
    when("A"   ){ $nh3 = $redn; }
    when("V"   ){ $voc = $redn; }
    when("S"   ){ $sox = $redn; }
    when("N"   ){ $nox = $redn; }
    when("P"   ){ $pm25 = $pmco = $redn; }
    when("NP"  ){ $nox = $pm25 = $pmco = $redn; }
    when("SNP" ){ $sox = $nox = $pm25 = $pmco =  $redn; }
    when("AN"  ){ $nh3 = $nox = $redn; }
    when("SNAV"){ $sox = $nox = $nh3 = $voc = $redn; }
    when("ALL" ){ $sox=$nox=$voc=$nh3=$testp=$co=$pm25=$pmco=$redn; }
   #when("BASE") then no change!
  }

  my $femisdat = slurp("$DataDir/femis.dat");
  die "ERROR!! No country Num for $cc!\n"
    unless defined(my $ccnum = $country_nums{$cc}) or defined($city_lonlat{$cc});

  # using 0 here as long as emissions are guaranteed to contain either
  # only anthropogenic or only natural emissions perl 'country'
  my $ss = 0; # 100 = antropogenic sectors (1..10)
              # 0 = all sectors
  $femisdat .= "$ccnum $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n"
    if defined($ccnum);
  my @cx_add=();
  given($cc){ # Add split countries
    when("DE" ){ push(@cx_add,( 9, 10)); }
    when("KZT"){ push(@cx_add,(53, 75)); }
    when("TMT"){ push(@cx_add,(77, 79)); }
    when("UZT"){ push(@cx_add,(76, 78)); }
    when("RU" ){ push(@cx_add,(36..38, 42)); }
    when("RUE"){ push(@cx_add,(36, 37, 38, 42, 71, 74)); } # Add RU and RFE
    when("ASI"){ push(@cx_add,(76, 77, 80, 82, 84)); }     # Add splitted ASI
    when("AST"){ push(@cx_add,(80, 82, 83, 84, 85)); }     # Add ASI and ASE
    when("ATL"){ push(@cx_add,(70)); }                     # Add ATL outside EMEP
              # Add split sea areas
    when("BAS"){ push(@cx_add,(302..309)); }
    when("NOS"){ push(@cx_add,(312..319)); }
    when("ATL"){ push(@cx_add,(322..329)); }
    when("MED"){ push(@cx_add,(332..339)); }
    when("BLS"){ push(@cx_add,(342..349)); }
  }           # Add citySR pairs
  @cx_add=(@{$city_lonlat{$cc}}) if defined($city_lonlat{$cc});
  
  foreach my $cx (@cx_add) {
    $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
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
