#!/usr/bin/perl

#Common script for Njord, Stallo and Titan.
#Choose $VILJE=1 or $STALLO=1 or $TITAN=1
#___________________________________________________________________
# queue commands for PBS

#Queue system commands start with #PBS (these are not comments!)
# lnodes= number of nodes, ppn=processor per node (max8 on stallo)
# Stallo, use ib for infiniband (fast interconnect).
# Ve/Vilje, use ppn=16, do not use :ib
#VePBS -lnodes=2:ppn=16
#Stallo
#PBS -lnodes=64
##PBS -lnodes=8:ppn=8:ib
# wall time limit of run
#PBS -lwalltime=07:30:00
# lpmeme=memory to reserve per processor (max 16GB per node)
#PBS -lpmem=1000MB
#make results readable for others:
#PBS -W umask=0022
# account for billing
#PBS -A nn2890k
# multiple tasks for paralel SR runs (one task per country)
##PBS -t 1-56
#___________________________________________________________________


#___________________________________________________________________
#Titan queue commands

#Queue system commands start with #SBATCH (these are not comments!)
#SBATCH --account=nn2890k
#SBATCH --nodes=8 --ntasks-per-node=8 --constraint=ib
#SBATCH --mem-per-cpu=2G --time=06:00:00
#SBATCH --job-name=emep
#SBATCH --output=job.%N.%j.out
#SBATCH  --error=job.%N.%j.err

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
#Tips:            SNYKOV or STALLO  NJORD            TITAN
#  submit job     qsub run.pl       llsubmit run.pl  sbatch run.sh
#  queue status   qstat -u $USER    llq              squeue -u $USER
#  job status                       checkjob 3456    checkjob 3456
#  kill job       qdel 3456         llcancel 3456    scancel 3456
#  submit multi-task SR job (one task per country)
#                 qsub -t 1-56 run.pl                arrayrun 1-56 run.sh
######################################################################

use 5.6.0;
use strict;
use warnings;
use File::Copy qw();
use File::Compare;

$| = 1; # autoflush STDOUT

#Choose one machine
my $VILJE=0;  #1 if Ve or Vilje is used
my $STALLO=1; #1 if stallo is used
my $TITAN=0;  #1 if titan is used

# -j4 parallel make with 4 threads
my @MAKE = ("gmake", "-j4", "MACHINE=snow");
   @MAKE = ( "make", "-j4", "MACHINE=vilje")  if $VILJE==1 ;
   @MAKE = ( "make", "-j4", "MACHINE=stallo") if $STALLO==1 ;
   @MAKE = ( "make", "-j4", "MACHINE=titan")  if $TITAN==1 ;
die "Must choose STALLO **or** VILJE **or** TITAN!\n"
  unless $STALLO+$VILJE+$TITAN==1;
my $MAKEMODE=0; #="EMEP2010";  # make EMEP2010, make SR-EMEP2010

my %BENCHMARK;
# OpenSource 2008
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2005,emis=>"Modrun07/OpenSourceEmis"           ,archive=>1,chem=>"EmChem03");
# Dave's preference for EMEP:
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2006,emis=>"Modrun10/EMEP_trend_2000-2008/2006",archive=>1,chem=>"EmChem09");
# EECCA Default:
   %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun11/EMEP_trend_2000-2009/2008",archive=>1,chem=>"EmChem09soa",make=>"EMEP");
# Status Runs:
#  %BENCHMARK = (grid=>"EECCA" ,year=>2007,emis=>"Modrun09/2009-Trend2007-CEIP") ;
#  %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun10/2010-Trend2008_CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2009,emis=>"Modrun11/2011-Trend2009-CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2010,emis=>"Modrun12/2012-Trend2010-CEIP",archive=>1,chem=>"EmChem09soa",make=>"EMEP2010");
# Alternative domains:
#  %BENCHMARK = (grid=>"TNO28" ,year=>2008,emis=>"emis_TNO28"         ,archive=>1);
#  %BENCHMARK = (grid=>"MACC02",year=>2008,emis=>"2008_emis_EMEP_MACC",archive=>1) ;
# Candidate for removal/update: Do not arrchive!
#  %BENCHMARK = (grid=>"EECCA" ,year=>2005,emis=>"Modrun11/EMEP_trend_2000-2009/2005");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2006,emis=>"Modrun11/EMEP_trend_2000-2009/2006");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2007,emis=>"Modrun11/EMEP_trend_2000-2009/2007");
if (%BENCHMARK) {
# $BENCHMARK{'archive'} = 1;                        # save summary info in $DataDir
  $BENCHMARK{'debug'} = $BENCHMARK{'archive'};      # chech if all debug flags are .false.
  $BENCHMARK{'ndx'}   = 8 if $BENCHMARK{'archive'}; # number of procesors in
  $BENCHMARK{'ndy'}   = 8 if $BENCHMARK{'archive'}; # x and y directions
# Default setting, if not previously specified
  $BENCHMARK{'chem'}  = "EmChem09soa"
    unless $BENCHMARK{'chem'};  # chemical mecanism, e.g. OpenSource 2008
  $BENCHMARK{'make'}  = ($BENCHMARK{'chem'} eq "EmChem09soa")?"EMEP":"all"
    unless $BENCHMARK{'make'};  # make target, e.g. Status 2010
  $MAKEMODE = $BENCHMARK{'make'} unless $MAKEMODE;
}

my $INERIS_FACS=0;  # Used for timefactors, and e,g TNOxx tests
my $SR= 0;     # Set to 1 if source-receptor calculation
               # check also variables in package EMEP::Sr below!!

my $CWF=0;     # Set to N for 'N'-day forecast mode (0 otherwise)
my ($CWFBASE, $CWFDAYS, $CWFMETV, @CWFDATE, @CWFDUMP, $eCWF) if $CWF;
if ($CWF) {
  chop($CWFBASE = `date +%Y%m%d`);   # Forecast base date     (default today)
       $CWFDAYS = $CWF;              # Forecast lenght indays (default $CWF)
       $CWFMETV = 0;                 # Met.UTC version        (default 'none')
       $CWFBASE = $ENV{"DATE"} if $ENV{"DATE"}; # Forecast base date, lenght
       $CWFDAYS = $ENV{"NDAY"} if $ENV{"NDAY"}; #  & MetUTC version can be passed
       $CWFMETV = $ENV{"UTC"}  if $ENV{"UTC"} ; #  as environment variables
       $CWFBASE = shift if @ARGV;    # Forecast base date, lenght
       $CWFDAYS = shift if @ARGV;    #  & MetUTC version can be passed
       $CWFMETV = shift if @ARGV;    #  as argument to script
  $CWF=($CWFMETV)?"$CWFMETV-$CWFBASE":"$CWFBASE";
  my $metday = ($CWFMETV)?int($CWFMETV/24+0.99):0;
  chop($CWFBASE = `date -d '$CWFBASE $metday day' +%Y%m%d`) if ($metday gt 0);
  chop($CWFDATE[0] = `date -d '$CWFBASE 1 day ago'    +%Y%m%d`);  # yesterday
       $CWFDATE[1] = $CWFBASE;                                    # start date
  chop($CWFDATE[2] = `date -d '$CWFBASE ($CWFDAYS-1) day' +%Y%m%d`);  # end date
  chop($CWFDUMP[0] = `date -d '$CWFBASE 1 day' +%Y%m%d000000`); # 1st dump/nest
  chop($CWFDUMP[1] = `date -d '$CWFBASE 2 day' +%Y%m%d000000`); # 2nd dump/nest
  $eCWF=0;                           # Emergency forecast
  $MAKEMODE=$eCWF?"eEMEP":"MACC";    # Standard Forecast model setup
# $MAKEMODE=$eCWF?"eEMEP2010":"MACC-EVA2010";    # 2010 special
}
 $CWF=0 if %BENCHMARK;
 $MAKEMODE="SR-$MAKEMODE" if($MAKEMODE and $SR);

# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables
#      are explained below, and derived variables set later.-

my $year = "2009";
   $year = substr($CWFBASE,0,4) if $CWF;
   $year = $BENCHMARK{"year"} if %BENCHMARK;
( my $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

my $iyr_trend = $year;
$iyr_trend = "2020" if $SR ;  # 2020 assumed for SR runs here
#$iyr_trend = "2020" ;  #  TSAP test

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
my $AGNES      = "nyiri";
my $ALVARO     = "alvarov";
my $ROBERT     = "mifarb";
my $HALDIS     = "mifahb";
my $BIRTHE     = "birthems";
my $FORCAST    = "forecast";

my $USER = $ENV{"USER"};
if ($PETER =~ m/$USER/) { $USER="$PETER" };
print "USER = $USER\n";


# hb NH3Emis
my $NH3EMIS_VAR = 0; # set to 1 if new temp NH3.

my $METformat="cdf"; # felt or cdf

my $GRID = "TNO28"; # TNO7, TNO14, TNO28, TNO56, EMEP, EECCA, MACC02, GLOBAL or FORECAST
   $GRID = $eCWF?"GLOBAL":"MACC02" if $CWF;
   $GRID = $BENCHMARK{'grid'} if %BENCHMARK;
#DS Confusing list of possibilites. Needs  CHECK LATER
my $MetDriver = "H20" ; # DS consider condition "EC";  #"H20";
   $MetDriver = "EC" if $year >= 2005; # Available Nov 2011

my ($HOMEROOT, $WORKROOT, $MetDir);
our $DataDir;
if ($STALLO) {
  $HOMEROOT = "/home";
  $WORKROOT = "/global/work";
  $DataDir  = "/global/work/$PETER/Data";
  $MetDir   = "$DataDir/$GRID/metdata/$year" ;
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year"  if ($GRID eq "GLOBAL");
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year"  if ($GRID =~ /TNO/);
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year"  if ($GRID eq "MACC02");
  $MetDir   = "$DataDir/$GRID/metdata_CWF/$year" if $CWF;
  $MetDir   = "$DataDir/$GRID/metdata_H20/$year" if $GRID eq "EECCA"; # assumes $METformat eq "cdf";

#DS added Mar 2011. Use EC for standard 2008 runs?
#Also print $MetDir into RunLog later
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year" if ($GRID eq "EECCA" && $year >= 2005 );

} elsif ($TITAN) {
  $HOMEROOT = "/usit/titan/u1";
  $WORKROOT = "/xanadu/d1";
  $DataDir  = "/projects/metno/emep/Data";
  $MetDir   = "$DataDir/$GRID/metdata/$year" ;
  $MetDir   = "$DataDir/$GRID/metdata_H20/$year" if $GRID eq "EECCA";
  $MetDir   = "$DataDir/$GRID/metcdf/$year" if ($GRID eq "EMEP") and
                                               ($METformat eq "cdf");
} else { #Ve or Vilje
  $HOMEROOT = "/home/metno";
  $WORKROOT = "$HOMEROOT/$USER/work";
  $DataDir  = "$HOMEROOT/mifapw/work/Data";
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year" ;
  $MetDir   = "$DataDir/$GRID/metdata_CWF/$year" if $CWF;
  $MetDir   = "/prod/forecast/work/emep/ec/prepmet" if $eCWF and ($USER eq $FORCAST);
}

# DataDir    = Main general Data directory
my $DATA_LOCAL = "$DataDir/$GRID";   # Grid specific data , EMEP, EECCA, GLOBAL
# Pollen data
my $PollenDir = "/home/$BIRTHE/Unify/MyData";
   $PollenDir = 0 unless $STALLO;
# Eruption (eEMEP)
my $EruptionDir = "$HOMEROOT/$ALVARO/Unify/MyData";
   $EruptionDir = 0 unless -d $EruptionDir;


# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

my $CityZen = 0 ;
  #$Chem     = "Eucaari_Trends";      # Label for chemical scheme used
my $VBS   = 0;
my $Chem     = "EmChem09soa";
#$Chem     = "CRI_v2_R5";
   $Chem     = $BENCHMARK{'chem'} if $BENCHMARK{'chem'};

my $testv = "rv4beta22";

#User directories
my $ProgDir  = "$HOMEROOT/$USER/Unify/Unimod.$testv";   # input of source-code
   $ProgDir  = "/prod/forecast/emep/eemep/src/Unimod.$testv" if $eCWF and ($USER eq $FORCAST);
my $ChemDir  = "$ProgDir/ZCM_$Chem";
my $Specials = "specials";  # default
#$Specials = "TSAP_Jul2012";  # used for TSAP runs in July 2012
# Check:
die "No ProgDir! $ProgDir\n" unless -d $ProgDir;
die "No ChemDir! $ChemDir\n" unless -d $ChemDir;


#---- emislist --------------------------------------------------------
open(EMIS,"<$ProgDir/CM_emislist.csv") or die "Need CM_emislist.cvs file!\n";
  my @emislist = split(/,/,<EMIS>);
  print "EMISLIST ", join(" ", @emislist ), "\n";
close(EMIS);
#----  chem packages  (e.g. EmChembase PMmass ) -----------------------
open(CHEM,"<$ProgDir/CM_chempackages.txt") or die "Need CM_emislist.cvs file!\n";
  my @packages = <CHEM> or die "Need CM_chempackage.txt!\n" ;
  print "CHEM packages:\n @packages\n";
close(CHEM);
#----------------------------------------------------------------------

# Check that the code directory has the chem files we want:
# Now use mk.GenChem, CM_ files erased from ZCM_ directories
#die "Mis-Match chemistry, Unimod.$testv Chem: $Chem" if
#  ( File::Compare::compare( "$ProgDir/CM_ChemSpecs_ml.f90" , "$ChemDir/CM_ChemSpecs_ml.f90"));

my $WORKDIR     = "$WORKROOT/$USER/$testv.$year";  # working and result directory
   $WORKDIR     = "$WORKROOT/$testv.$year" if($WORKROOT =~ /$USER/);  # working and result directory
   $WORKDIR =~ s/$testv.$year/Benchmark\/$GRID.$year/g if (%BENCHMARK);
   $WORKDIR = "/prod/forecast/run/eemep" if $eCWF and ($USER eq $FORCAST);
my $MyDataDir   = "$HOMEROOT/$USER/Unify/MyData";           # for each user's private input
my $SoilDir     = "$DATA_LOCAL/dust_input";               # Saharan BIC
$SoilDir = 0 unless -d "$SoilDir/BC_DUST/2000";

# TEST! Road dust NOTE! The road dust code may not be working properly yet! Not tested enough!
my $RoadDir     = "$HOMEROOT/$ROBERT/Unify/MyData/TNO_traffic/" ;
$RoadDir = 0 if $CWF;
$RoadDir = 0 unless -d $RoadDir;#default?
#$RoadDir = 0 ;#default?


# Forecast: nest/dump dir, BCs pattern
my ($CWFDUMPDIR, $CWFBC) if $CWF;
if ($CWF) {
# $CWFDUMPDIR = "$WORKROOT/$USER/$testv.dump";
  ($CWFDUMPDIR = $WORKDIR) =~ s/$testv.$year/$testv.dump/g;
  $CWFBC = "$DataDir/$GRID/Boundary_conditions/%04d_IFS-MOZART_FC/cwf-mozifs_h%08d00_raqbc.nc" # IFS-MOZ Forecast
# $CWFBC = "$DataDir/$GRID/Boundary_conditions/%04d_IFS-MOZART_AN/h%08d00_raqbc.nc";           # IFS-MOZ ReAnalysus
# $CWFBC = "$DataDir/$GRID/Boundary_conditions/%04d_EVA/EVA_%08d_EU_AQ.nc" # IFS-MOZ Forecast
}

#ds check: and change
chdir "$ProgDir";
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"
#  unless $ENV{PWD} =~ /Unimod.$testv.$year/;
print "TESTING ENV:", $ENV{PWD}, "\n";


# Default emissplits used here. if $Specials is set will look
my $SplitDir    = "$DataDir/SPLITS_JAN2010/BASE_NAEI2000_GH2009.$Chem" ;
$SplitDir    = "$ChemDir/EMISSPLIT";
#RB:had "~mifarb/Unify/MyData/D_EGU/SPLITS_NOV2009v2/BASE_NAEI2000_GH2009.$Chem" ;

my $version     = "Unimod" ;
my $PROGRAM     = "$ProgDir/$version";         # programme
my $subv        = "$testv" ;                  # sub-version (to track changes)

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later.
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR

#Possible emission scenarios for HIRHAM run
#GEA scenarios: HIGH_CLE, HIGH_FROZEN, LOW_CLE, LOW_SLE
#Historic emissions from Lamargue et al. : historic_emis
my $emisscen = "historic_emis";
my $emisyear = $year;
$scenario = "${emisscen}_emis${emisyear}_met${year}" if $GRID eq "HIRHAM";

my @runs     = ( $scenario );

#EMISSIONS: default settings
my ($EMIS_INP, $emisdir, $pm_emisdir);
$EMIS_INP = "$DATA_LOCAL"                            ;
$EMIS_INP = "$DATA_LOCAL/Emissions/Modruns" if $TITAN;

#dave: Use Modrun11 if possible:
my $EMIS_OLD = "/global/work/$AGNES/Emission_Trends";
$emisdir = "$EMIS_OLD/$year" if $year < 2000;
$emisdir = "$EMIS_INP/Modrun11/EMEP_trend_2000-2009/$year" if ( $year > 1999 ) and ($year < 2009);
$emisdir = "$EMIS_INP/Modrun11/2011-Trend2009-CEIP" if $year >= 2009 ;

#Alvaro: When are we to switch to Modrun12?
#$emisdir = "$EMIS_INP/Modrun12/EMEP_trend_2000-2009/$year" if ( $year > 1999 ) and ($year < 2010);
#$emisdir = "$EMIS_INP/Modrun12/2012-Trend2010-CEIP" if $year >= 2010 ;

#TMP and should be improved because it gives errors for
# other domains!
#.. For using emissions of EC/OC instead of PMx
my $RFEmisDir = "/global/work/$SVETLANA/Data_RF"; # Split-Fraction files for EC/OC
my $TNOemisDir = "/global/work/$SVETLANA/Emis_TNO"; # TNO EC/OC emissions

#$emisdir = $TNOemisDir if $EUCAARI;
$emisdir = "$EMIS_INP/emissions/${emisscen}/${emisyear}" if $GRID eq "HIRHAM";

$pm_emisdir = $emisdir;
$pm_emisdir = "$EMIS_INP/2006-Trend2000-V7"  if $year < 2000;
$pm_emisdir = "/home/$ROBERT/Unify/MyData/D_EGU/${GRID}_GRID" if $VBS;

#EMISSIONS: FORECAST settings
if ( ($GRID eq "FORECAST") or ($GRID eq "GEMS025") or ($GRID eq "MACC02") ) {
  $EMIS_INP = "$DATA_LOCAL/Emissions";
  $emisdir = "$EMIS_INP/2008-Trend2006-V9-Extended_PM_corrected-V3"; # GEMS025
# $emisdir = "$EMIS_INP/2008_emis_EMEP_from_PS50" if $GRID eq "MACC02"; # MACC02
# $emisdir = "$EMIS_INP/2008_emis_EMEP_MACC" if $GRID eq "MACC02"; # MACC02
  $emisdir = "$EMIS_INP/2007_emis_MACC" if $GRID eq "MACC02"; # MACC02
  $pm_emisdir = $emisdir;
}
if ( $GRID =~ /TNO/)  {  # STALLO ONLY
  die "Need to set emisdir for TNO, non-stallo" unless $STALLO;
  $EMIS_INP = "/global/work/nyiri/emis_SRbase/INERIS_direct/$GRID";
  $emisdir = $EMIS_INP;
  $pm_emisdir = $emisdir;
}

#EMISSIONS: BENCHMARK settings
if (%BENCHMARK){
  $emisdir    = "$EMIS_INP/$BENCHMARK{'emis'}";
  $pm_emisdir = $emisdir;
}

#Dave, reset to Emission_Trends for Chem project, Oct 18th
my $TREND_RUNS = 0;
if ($STALLO && $TREND_RUNS ) {
  $EMIS_INP = "/global/work/$AGNES/Emission_Trends";
  die "Year not in trend run series!! " unless -f $EMIS_INP/$year;
  $emisdir = "$EMIS_INP/$year";
  $pm_emisdir = $emisdir;
}
#
if ( $GRID eq "GLOBAL" ) {
  $EMIS_INP = $DATA_LOCAL;
  $emisdir = ($eCWF)?"$EMIS_INP/Emissions_June2012":"$EMIS_INP/MonthlyEmis";
  $pm_emisdir = $emisdir;
}

my $RESET        = 0 ;  # usually 0 (false) is ok, but set to 1 for full restart
my $COMPILE_ONLY = 0 ;  # usually 0 (false) is ok, but set to 1 for compile-only
my $INTERACTIVE  = 0 ;  # usually 0 (false), but set to 1 to make program stop
my $DRY_RUN      = 0 ;  # Test script without running model (but compiling)

# just before execution - so code can be run interactivel.

# NDX, NDY  now set in ModelConstants_ml - use perl to extract these
# values and check against submission:

print "$ProgDir/ModelConstants_ml.f90\n";
open(IN,"<$ProgDir/ModelConstants_ml.f90");
my ( $NDX, $NDY , $XDIM ); # Processors in x-, y-, direction, and x-dimension
while(my $line = <IN>){
    $line=~ s/!.*//; # Get rid of comment lines
    $NDX = $1 if $line =~ /\W+ NPROCX \W+ (\d+) /x ;
    $NDY = $1 if $line =~ /\W+ NPROCY \W+ (\d+) /x ;
    $XDIM = $1 if $line =~ /\W+ IIFULLDOM \W+ (\d+) /x ;
}
close(IN);
my $NPROC =  $NDX * $NDY ;
print "ModelConstants has: NDX = $NDX NDY = $NDY  =>  NPROC = $NPROC\n";
die "Global model requires NDY <= 2\n" if ( $GRID eq "GLOBAL" && $NDY > 2);
die "Domain mis-match Model: $XDIM Grid $GRID" if (
   ( $GRID eq "EMEP"     && $XDIM != 170 ) or
   ( $GRID eq "EECCA"    && $XDIM != 132 ) or
   ( $GRID eq "TNO7"     && $XDIM != 840 ) or
   ( $GRID eq "TNO14"    && $XDIM != 420 ) or
   ( $GRID eq "TNO28"    && $XDIM != 210 ) or
   ( $GRID eq "TNO56"    && $XDIM != 105 ) or
   ( $GRID eq "EECCA_25" && $XDIM != 264 ) or
   ( $GRID eq "MACC02"   && $XDIM != 321 ) or
   ( $GRID eq "HIRHAM"   && $XDIM != 182 ) or
   ( $GRID eq "GLOBAL"   && $XDIM != 360 ) );

if (%BENCHMARK and $BENCHMARK{'debug'}){
  die "No debug flags for benchmarks!"
  if (system ("grep -Hnie 'logical.*DEBUG.*=\ *.TRUE.' $ProgDir/*.f90") == 0) or
     (system ("grep -Hnie 'DEBUG.*=\ *.TRUE.' $ProgDir/ModelConstants_ml.f90") == 0);
}

if ( $ENV{PBS_NODEFILE} ) {
   $_ =  `wc -l $ENV{PBS_NODEFILE}`;
   my $RUN_NPROC;
   ( $RUN_NPROC, undef ) = split;
   print "Qsub has: lnodes $RUN_NPROC\n";
   die "Error: Wrong number of lnodes!\n" unless $NPROC == $RUN_NPROC;
   die "Domain mis-match Model: ${NDX}x${NDY} for Benchmark"
   if(($BENCHMARK{'ndx'} and $BENCHMARK{'ndx'} != $NDX) or
      ($BENCHMARK{'ndy'} and $BENCHMARK{'ndy'} != $NDY));
} else {
   print "skip nodefile check on interactive runs\n";
}

# P B S -lnodes=$NPROC
#(Didn't work: my $NCPUS=`wc -l $ENV{PBS_NODEFILE} | awk '{print $1}'`;)



my @month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

#Only 360 days in HIRHAM metdata. We ignore leaps
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,24) if $GRID eq "HIRHAM";

my $mm1   =  "07";       # first month, use 2-digits!
my $mm2   =  "07";       # last month, use 2-digits!
my $dd1   =  1;       # Start day, usually 1
my $dd2   =  1;       # End day (can be too large; will be limited to max number of days in the month)
                      # put dd2=0 for 3 hours run/test.

if (%BENCHMARK){ # Allways runn full year on benchmark mode
  $mm1   =  "01";
  $mm2   =  "12";
  $dd1   =  1;       # Start day, usually 1
  $dd2   =  31;       # End day, usually 31
}

# <---------- end of normal use section ---------------------->
# <---------- end of user-changeable section ----------------->
#               (normally, that is...)

if ($SR) {
    print "SR is true\n";
    @runs = EMEP::Sr::initRuns();
}

if ($CWF) {
  print "CWF is true: $CWFDAYS-day foracast mode\n";
  @runs = ( $CWF ) unless $SR ;
}

if (%BENCHMARK){
  print "BENCHMARK is true: $GRID $year $testv $Chem\n";
  @runs = ( "$testv-$Chem" );
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
  if ($MAKEMODE) {
    system(@MAKE, "$MAKEMODE");
  } elsif ($SR) {
    #No recompile in SR runs
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

if ($MAKEMODE) {
  system(@MAKE, "$MAKEMODE") == 0 or die "@MAKE $MAKEMODE failed";
} else {
  system (@MAKE, "depend") ;
  system (@MAKE, "all") == 0 or die "@MAKE all failed";
}

die "Done. COMPILE ONLY\n" if  $COMPILE_ONLY;  ## exit after make ##


my @list_of_files = ();   # Keep list of data-files
my $cwfbc = "No BC file";


########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################
########################### START OF RUNS  ##########################

foreach my $scenflag ( @runs ) {
  if ($SR) {
    $scenario = EMEP::Sr::getScenario($scenflag);
  } elsif ($CWF) {
    $scenario = $eCWF?"eemep-$scenflag":"CWF_$scenflag";
  } elsif (%BENCHMARK){
    $scenario = "BM_$scenflag";
  } else {
    $scenario = $scenflag;
  }
  print "STARTING RUN $scenario \n";

  my $runlabel1    = "$scenario";   # NO SPACES! SHORT name (used in CDF names)
  my $runlabel2    = "$testv\_$Chem\_$scenario\_$year\_Trend$iyr_trend";   # NO SPACES! LONG (written into CDF files)

  my $RESDIR = "$WORKDIR/$scenario";
  mkdir_p($RESDIR);

  chdir $RESDIR;   ############ ------ Change to RESDIR
  print "Working in directory: $RESDIR\n";

  if ($METformat eq "felt") {
    my $nnn = 1;
    for (my $mm = $mm1; $mm <= $mm2; $mm++) {
      # Assign met files to fil001...etc.
      for (my $n = 1; $n <= $month_days[$mm]; $n++) {
        $nnn = metlink($n, $nnn, $mm);
      }
    }

    my $mmlast = $mm2 + 1;
    my $yylast = $year;
   if ( $mmlast > 12 ) {
      $yylast = $yylast + 1;
      $mmlast = 1;
    }
    my $old = sprintf "$MetDir/f00.%04d%02d01", $yylast, $mmlast;
    my $new = sprintf "fil%04d", $nnn;
    mylink( "LAST RECORD SET: ", $old,$new ) ;

  } elsif ($CWF) {
# Forecast Meteorology in NetCDF
    for (my $n = 0; $n < $CWFDAYS; $n++) {
      my $metday = int($CWFMETV/24+0.99);
     (my $old = "$MetDir/meteo%Y%m%d${CWFMETV}_%%02d.nc") =~ s/$year/%Y/g;
      chop($old=`date -d '$CWFBASE $metday day ago' '+$old'`);
      $old = sprintf "$old",$n+$metday;
      if (-e $old) {
        chop($CWFDATE[2]=`date -d '$CWFBASE $n day' +%Y%m%d`);
        my $new = "meteo$CWFDATE[2].nc";
        mylink( "Linking:", $old,$new ) ;
      # IFS-MOZART BC file
        $old = sprintf "$CWFBC",substr($CWFDATE[2],0,4),$CWFDATE[2];
        $new = "EMEP_IN_BC_$CWFDATE[2].nc";
        if (-e $old) {
          mylink( "Linking:", $old,$new );
          $cwfbc=$old;
        } else {
          die "BC file for $CWFDATE[2] not available (yet)" if ($CWFBC =~ m/_AN.+h$year/);
          print "BC file for $CWFDATE[2] not available (yet). Try $CWFDATE[0] BC file\n";
          $old = sprintf "$CWFBC",substr($CWFDATE[0],0,4),$CWFDATE[0];
          if (-e $old) {
            mylink( "Linking:", $old,$new ) if (-e $old);
            $cwfbc=$old;
          }
        }
      } else {
      # meteo not in place !!!!!
        die "Meteo file $old for $CWFBASE not available (yet). Try later ...\n";
      }
    }
# Forecast nest/dump files
    my $old="$CWFDUMPDIR/CWF_${CWFDATE[0]}_dump.nc";  # yesterday's BASE dump
      ($old="$CWFDUMPDIR/${scenario}_dump.nc")        # today's dump
            =~s/$CWFBASE/$CWFDATE[0]/g;               # yesterday's dump
    if (-e $old) {
      # we manage to link the dumpfile
      my $new="EMEP_IN_IC.nc";
      mylink( "Linking:", $old, $new);
    } else {
      # if we have yesterday meteo, we can have an extra spin up day!
      print "No dumpfile present:\n\t$old,\n\ttrying to take an extra spin-up day!\n";
      my $old = "$MetDir/meteo${CWFDATE[0]}_00.nc";
      if (-e $old) {
        my $new = "meteo${CWFDATE[0]}.nc";
        mylink( "Linking:", $old, $new);
        print "Managed to link meteo for an extra spin-up day!\n";
      # Update simulation: start-date ($CWFDATE[1]), "yesterday" ($CWFDATE[0])
        $CWFDATE[1]=$CWFDATE[0];
        chop($CWFDATE[0]=`date -d '$CWFDATE[0] 1 day ago' +%Y%m%d`);
      # IFS-MOZART BC file
        $old = sprintf "$CWFBC",substr($CWFDATE[1],0,4),$CWFDATE[1];
        $new = "EMEP_IN_BC_${CWFDATE[1]}.nc";
        if (-e $old) {
      # we manage to link the BC file
          mylink( "Linking:", $old,$new );
          print "Managed to link BC file for the extra spin-up day!\n";
          $cwfbc=$old;
        } else {
          print "BC file for $CWFDATE[1] not available (yet) (spin-up). Try yesterdays BC file\n";
          $old = sprintf "$CWFBC",substr($CWFDATE[0],0,4),$CWFDATE[0];
          if (-e $old) {
            mylink( "Linking:", $old,$new );
            print "Managed to link yesterdays BC file for the extra spin-up day!\n";
            $cwfbc=$old;
          }
        }
      # see if we can link to a dump file ...
        my $old="$CWFDUMPDIR/CWF_${CWFDATE[0]}_dump.nc"; # yesterday's BASE dump
          ($old="$CWFDUMPDIR/${scenario}_dump.nc")       # today's dump
                =~s/$CWFBASE/$CWFDATE[0]/g;              # yesterday's dump
        if (-e $old) {
        # we manage to link the dumpfile
          my $new="EMEP_IN_IC.nc";
          mylink( "Linking:", $old, $new);
          print "Managed to link dump file for the extra spin-up day!\n";
        }
      }
    }
# Update start and end months (used for linking some climatological files)
    $mm1=substr($CWFDATE[1],4,2);  # start date
    $mm2=substr($CWFDATE[2],4,2);  # end date
  } else {
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
    if ( $mmlast > 12 ) {
      $yylast = $yylast + 1;
      $mmlast = 1;
    }
    my $old = sprintf "$MetDir/meteo%02d%02d01.nc", $yylast, $mmlast;
    my $new = sprintf "meteo%02d%02d01.nc", $yylast, $mmlast;
    mylink( "LAST RECORD SET: ", $old,$new ) ;
  }

#=================== INPUT FILES =========================================
# ToDo Change noxsplit.default to defaults, as with voc (also in Unimod)

  my %ifile   = ();   # List of input data-files

# First, emission files are labelled e.g. gridSOx, which we assign to
# emislist.sox to ensure compatability with the names (sox,...) used
# in the model. It doesn't matter if we have extra mapping here,it
# is  GenIn.reactions and the associated emislist that decides what gets used.
# e.g. lines such as:
#  emisfiles:sox,nox,co,voc,nh3
#  emisfiles:pm25
# etc.

my $timeseries  = "$DataDir/inputs_emepdefaults_Jun2012";
my $Tbase = 18 ;  # Base-temperature for Degree-day (HDD) files

if ( $INERIS_FACS  ){
   $timeseries  = "$DataDir/inputs_eurodelta_Jun2012";
   $Tbase = 20;
}
   #$ifile{"$timeseries/HourlyFacs.EMEP2003"} = "HOURLY-FACS";
   #$ifile{"$timeseries/HourlyFacs.TNO2005"} = "HOURLY-FACS";
# INERIS provided the most complete hourly file, we use as default
$ifile{"$timeseries/HourlyFacs.INERIS"} = "HOURLY-FACS";

$ifile{"$timeseries/EmisHeights.txt"} = "EmisHeights.txt";

  my %gridmap = ( "co" => "CO", "nh3" => "NH3", "voc" => "NMVOC",
                  "sox" => "SOx", "nox" => "NOx" ,
                  "pm10" => "PM10", "pm25" => "PM25", "pmco" => "PMco",
 # VBS specials
                  "pocfwd" => "POCfWD",
                  "pocffl" => "POCfFL", "poccfl"   => "POCcFL",
                  "ecfwd" => "ECfWD", "eccwd" => "ECcWD",
                  "ecffl" => "ECfFL", "eccfl" => "ECcFL",
                  "forfbc"   => "FORFBC", "forfoc"   => "FORFOC",
 #  Sometimes used also:
                  "ecfi" => "ECfine","ecco" => "ECcoar", "ocfi" => "OCfine" ) ;
                  # sometimes was "ocfi" => "POCfine"   ) ;


  foreach my $poll  ( @emislist  ) {
    my $dir = $emisdir;
    $dir = $pm_emisdir if $poll =~ /pm/;   # FIX needed prior to 2000
 # VBS specials #rb Wood burning, Fossil fuel and Forest fire PM from TNO files
    $dir = $pm_emisdir if $poll =~ /wd/;   #
    $dir = $pm_emisdir if $poll =~ /fl/;   #
    $dir = $pm_emisdir if $poll =~ /forf/;   #
print "TESTING PM $poll $dir\n";

# hb NH3emis, new emis files
    if(($NH3EMIS_VAR)&&($poll eq "nh3")){
      $dir = "/home/$AGNES/emis_NMR";
      $ifile{"$dir/gridNH3_NMR_$year"} = "emislist.$poll";
    }else{
      $ifile{"$dir/grid$gridmap{$poll}"} = "emislist.$poll"
    }

    # copy pm25 if needed, avoid having 20 different PM25 time-series

      if ( -f "$timeseries/MonthlyFac.$poll" ) {
         print "FINDS??? Daily Fac pm25 fill in for $poll\n";
         system("wc $timeseries/MonthlyFac.$poll");
         $ifile{"$timeseries/MonthlyFac.$poll"} = "MonthlyFac.$poll";
         $ifile{"$timeseries/DailyFac.$poll"} = "DailyFac.$poll";
      } else { # Assume same as PM25, works e.g for ocffl, etc.
         print "Monthly Daily Fac pm25 fill in for $poll\n";
         my $TMPWDIR = "$WORKDIR/$testv.tmpdir";
         mkdir($TMPWDIR) unless -d $TMPWDIR;
         cp ("$timeseries/MonthlyFac.pm25", "$TMPWDIR/MonthlyFac.$poll");
         cp ("$timeseries/DailyFac.pm25",   "$TMPWDIR/DailyFac.$poll");
         $ifile{"$TMPWDIR/MonthlyFac.$poll"} = "MonthlyFac.$poll";
         $ifile{"$TMPWDIR/DailyFac.$poll"} = "DailyFac.$poll";

      }

    $ifile{"$SplitDir/emissplit.defaults.$poll"} = "emissplit.defaults.$poll";
    # specials aren't essential, but if available we use them
    # Set $Specials flag for special cases, e.g. TSAP
    # INERIS special! nox and pm. Take from 2010 IIASA
    #if ( $INERIS_FACS && -e "$timeseries/emissplit.specials.$poll.2010" ) {

    if ( -e "$timeseries/emissplit.$Specials.$poll.$iyr_trend" ) { # e.g. TSAP
        $ifile{"$timeseries/emissplit.$Specials.$poll.$iyr_trend"} =
               "emissplit.specials.$poll"
    } elsif ( -e "$timeseries/emissplit.specials.$poll.2010" ) { # when no other year availanle
        $ifile{"$timeseries/emissplit.specials.$poll.2010"} =
               "emissplit.specials.$poll"
    } else {

     $ifile{"$SplitDir/emissplit.specials.$poll"} = "emissplit.specials.$poll"
     if( -e "$SplitDir/emissplit.specials.$poll" );
    }
  }

  foreach my $mmm ( $mm1 .. $mm2, $mm1, $mm2 ) {
    my $mm = sprintf "%2.2d", $mmm;
    $ifile{"$DATA_LOCAL/natso2$mm.dat"} =  "natso2$mm.dat";
    $ifile{"$DataDir/lt21-nox.dat$mm"} =  "lightning$mm.dat";
# BIC for Saharan dust
    if ( $SoilDir ) { # Not yet for EMEP domain
      foreach my $bc ( qw ( DUST_c_ext DUST_f_ext )) { #
        $ifile{"$SoilDir/BC_DUST/2000/$bc.$mm"} =  "$bc.$mm";
      }
    } # dust
    if ( $GRID eq "GLOBAL" ) {
      $mm = sprintf "%2.2d", $mmm;
      foreach my $t ( qw (nox voc co nh3 pm25 pmco) ) {
        $ifile{"$emisdir/grid$gridmap{$t}.$mm"} =  "grid$t.$mm";
      }
      foreach my $t ( qw (so2) ) {
        $ifile{"$emisdir/gridSO2.$mm"} =  "gridsox.$mm";
      }
    }
  }

# Emissions setup:
#  if ($EUCAARI) { # DS RE-CHECK shouldn't be needed
#    $ifile{"$TNOemisDir/femis.dat"} =  "femis.dat";
#    $ifile{"$DATA_LOCAL/emissions/femis.dat"} =  "femis.dat" if $GRID eq "HIRHAM" ;
#  } else {
    $ifile{"$ChemDir/femis.defaults"} =  "femis.defaults";  # created now by GenChem

# my $old="$DATA_LOCAL/Boundary_and_Initial_Conditions.nc";
# my $new="Boundary_and_Initial_Conditions.nc";
# mylink( "BIC: ", $old,$new ) ;
#EUCAARI, but all?
# Skip:  $ifile{"$DATA_LOCAL/Boundary_and_Initial_Conditions.nc"} =
#                     "Boundary_and_Initial_Conditions.nc" unless $GRID eq "MACC02";
  $ifile{"$DataDir/GLOBAL_O3.nc"} =
                  "GLOBAL_O3.nc";
  $ifile{"$DataDir/amilt42-nox.dat"} = "ancatmil.dat";#RENAME TO AIRCARAFT?!
  $ifile{"$DataDir/GLOBAL_ForestFireEmis.nc"} = "GLOBAL_ForestFireEmis.nc"; #GFED emissions
  $ifile{"$DataDir/ForestFire_Emis_$year.nc"} = "GLOBAL_ForestFireEmis_FINN.nc"
    if ($year >= 2002 and $year <= 2011);#FINN emissions
  $ifile{"$DataDir/nox_emission_1996-2005.nc"} = "nox_emission_1996-2005.nc";
  $ifile{"$DataDir/AircraftEmis_FL.nc"} = "AircraftEmis_FL.nc";
  $ifile{"$DataDir/SurfacePressure.nc"} = "SurfacePressure.nc";
  $ifile{"$DataDir/SoilTypes_IFS.nc"} = "SoilTypes_IFS.nc";
#netcdf RoadDust inputs:
  $ifile{"$DataDir/RoadMap.nc"} = "RoadMap.nc";
  $ifile{"$DataDir/AVG_SMI_2005_2010.nc"} = "AVG_SMI_2005_2010.nc";


#TEMPORARY SETUP
#  my $tmpndep = "/home/$DAVE/Work/RESULTS/MAPS/AnnualSums/AnnualNdep";
#  $ifile{"$tmpndep/AnnualNdep_BM_rv3_9_20soa-EmChem09soa.nc"} = "AnnualNdep.nc";

if ( $iyr_trend > 2015 )  {
  $ifile{"$DataDir/AnnualNdep_TNO28_2020.nc"} = "annualNdep.nc";
} else { # Note PS50x - hand-edited version
  $ifile{"$DataDir/AnnualNdep_PS50x_EECCA2005_2009.nc"} = "annualNdep.nc";
}

# hb NH3emis
# New ammonia emissions  ---   NB no read permissions yet!!
  $ifile{"/home/$HALDIS/Unimod_NMR_NH3/Unimod.rv3_6_8/Sector_NH3Emis.txt"}="Sector_NH3Emis.txt" if($NH3EMIS_VAR);

# new inputs style (Aug 2007)  with compulsory headers:
# From rv3_14 used only for FORECAST mode
  $ifile{"$DATA_LOCAL/Inputs.Landuse"} = "Inputs.Landuse" if ( $CWF ) ;
  $ifile{"$DataDir/Landuse/landuseGLC2000_INT1.nc"} ="GLOBAL_landuse.nc";

  $ifile{"$DataDir/Landuse_PS_5km.nc"} ="Landuse_PS_5km.nc";
  $ifile{"$DataDir/LanduseGLC.nc"} ="LanduseGLC.nc";

  $ifile{"$DataDir/LandInputs_Jan2012/Inputs_DO3SE.csv"} = "Inputs_DO3SE.csv";
  $ifile{"$DataDir/LandInputs_Jan2012/Inputs_LandDefs.csv"} = "Inputs_LandDefs.csv";

  $ifile{"$DataDir/sondesLL.dat"} = "sondes.dat";
  $ifile{"$DataDir/sitesLL.dat"} = "sites.dat";
  # Extended to get isoprene, HCHO EC, OC /(huge list!)
  # not default
  #$ifile{"$MyDataDir/sitesCPM_ds.dat"} = "sites.dat";

# DEGREE DAYS (Tbase set above, either 18 or 20):
#
 my $HDD = "$timeseries/HDD${Tbase}-${GRID}-$year.nc";
 print "Looking for DegreeDayFac: $HDD \n";
 my $skipHDD = 0;
 $skipHDD = 1 if $CWF ;
 $skipHDD = 1 if $GRID eq "RCA" ;
 $skipHDD = 1 if $GRID eq "GLOBAL" ;
 unless ( $skipHDD ) {
   system("ls -lh --time-style=long-iso $HDD");
   die "NO HDD files " unless -f $HDD;   # Can comment out if USE_DEGREEDAYS
                                      # set false in ModelConstants_ml
   $ifile{"$HDD"} = "DegreeDayFactors.nc" if -f $HDD ;
 }

  $ifile{"$DataDir/GLOBAL_LAInBVOC.nc"} = "GLOBAL_LAInBVOC.nc";
#New EURO BVOC
  $ifile{"$DataDir/LandInputs_Mar2011/EMEP_EuroBVOC.nc"} = "EMEP_EuroBVOC.nc";

# Seasonal stuff  ----    Can't we improve this? e.g. every month?
  my %seasons = ( "jan" => "01", "apr" => "02", "jul" => "03" , "oct"=> "04");
  foreach my $s ( keys(%seasons) ) {
    $ifile{"$DataDir/a${s}t42-nox.dat"} = "ancat$seasons{$s}.dat";
    $ifile{"$DataDir/jclear.$s"} = "jclear$seasons{$s}.dat";
    $ifile{"$DataDir/jcl1.$s"} = "jcl1km$seasons{$s}.dat";
    $ifile{"$DataDir/jcl3.$s"} = "jcl3km$seasons{$s}.dat";
  }

#  $ifile{"$DATA_LOCAL/rough.dat"} = "landsea_mask.dat"; # Roughness length;
  #NOTNEEDED $ifile{"$DATA_LOCAL/Volcanoes.dat"} = "Volcanoes.dat" unless $EUCAARI;
  $ifile{"$DataDir/VolcanoesLL.dat"} = "VolcanoesLL.dat";
# Volcanic Eruption (eEMEP)
  if(($MAKEMODE =~ /2010/) or ($MAKEMODE =~ /eEMEP/)){
    cp ("$ChemDir/eruptions.csv","eruptions.csv");
    $ifile{"$ChemDir/volcanoes.csv"} = "volcanoes.csv";
    print "$ChemDir/volcanoes.csv\n";
    open(IN,"<$ChemDir/volcanoes.csv");
    while(my $line = <IN>){
      $line=~ s/#.*//;                 # Get rid of comment lines
      my $vname = (split(",",$line))[0];  # Volcano tracer name
      my $efile = "$EruptionDir/${vname}_7bin.eruptions";
      $efile = "$EruptionDir/${vname}_2bin_${MAKEMODE}.eruptions" if($MAKEMODE =~ /2010/);
    # $efile = "$EruptionDir/${vname}_2bin_${MAKEMODE}_SR-SOx.eruptions" if($MAKEMODE =~ /2010/);
    # $efile = "$EruptionDir/${vname}_2bin_${MAKEMODE}_SR-PMx.eruptions" if($MAKEMODE =~ /2010/);
      system("cat $efile >> eruptions.csv") if -e $efile;
    }
    close(IN);
  }

# For Pollen
  if ( $PollenDir ) {
      $ifile{"$PollenDir/pollen_data.nc"} = "pollen_data.nc";
  }


# For windblown dust
   if ( $SoilDir ) {
    $ifile{"$SoilDir/clay_isric_percent_ext.dat"} = "clay_frac.dat";
    $ifile{"$SoilDir/sand_isric_percent_ext.dat"} = "sand_frac.dat";
   }

# TEST!!! Road dust NOTE! The road dust code is NOT thoroughly tested yet!
# NOTE ALSO THAT the Climate factors in the file below are just rough estimates based on the TNO soil water data, to be updated with something based on EMEP soil water!
  if ( ($GRID eq "EECCA") and ($RoadDir) ) {
      $ifile{"$RoadDir/RoadDust_HIGHWAYplus_emis_potential.txt"} = "HIGHWAYplus";
      $ifile{"$RoadDir/RoadDust_NonHighway_emis_potential.txt"} = "NONHIGHWAY";
      $ifile{"$RoadDir/RoughTestClimateFactorSoilWater.txt"} = "ROADDUST_CLIMATE_FAC";
  }elsif( $GRID =~ /TNO/ ){
      $ifile{"$DATA_LOCAL/nonHIGHWAYs_RoadDust_potentials.txt"} = "NONHIGHWAY";
      $ifile{"$DATA_LOCAL/HIGHWAYplus_RoadDust_potentials.txt"} = "HIGHWAYplus";
      $ifile{"$DATA_LOCAL/ClimateFactors_SMI.txt"} = "ROADDUST_CLIMATE_FAC";
  }
# IFZ-MOZ BCs levels description (in cdo zaxisdes/eta format)
  $ifile{"$DataDir/$GRID/Boundary_conditions/mozart_eta.zaxis"} = "EMEP_IN_BC_eta.zaxis"
    if ( $CWF and -e $cwfbc and $cwfbc =~ m/IFS-MOZART/ );

  foreach my $old ( sort keys %ifile ) {  # CHECK and LINK
    if ( -r $old ) {
      my $new =  $ifile{$old};
      mylink( "Inputs: ", $old,$new ) ;
    } else {
      print "Missing Input $old !!!\n";
      die "ERROR: Missing OLD $old (or possibly wrong Chem$Chem\n" unless $old =~ /special/;
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

# Link executable also, since gridur is funny about these
# things....

  my $LPROG = "Unimod";
 #mylink( "PROGRAM!!  ", $PROGRAM,$LPROG) ;
  cp ($PROGRAM, $LPROG) or die "cannot copy $PROGRAM to $LPROG\n";
  push(@list_of_files , $LPROG);    # For later deletion

# Write out list of linked files to a shell-script, useful in case the program
# hangs or crashes:

  open(RMF,">Remove.sh");
  foreach my $f ( @list_of_files ) { print RMF "rm $f \n"; }
  print RMF "rm $LPROG\n";   # Delete executable also
  close(RMF);

  my $startdate = sprintf "%04d %02d %02d",$year,$mm1,$dd1;
  my $enddate   = sprintf "%04d %02d %02d",$year,$mm2,$dd2;
  if ($CWF){    # use forecast start date
    $startdate = substr($CWFDATE[1],0,4)." ".substr($CWFDATE[1],4,2)." ".substr($CWFDATE[1],6,2);
    $enddate   = substr($CWFDATE[2],0,4)." ".substr($CWFDATE[2],4,2)." ".substr($CWFDATE[2],6,2);
  }

# make file with input parameters (to be read by Unimod.f90)
  unlink("INPUT.PARA");
  open(TMP,">INPUT.PARA");
  print TMP "$iyr_trend\n$runlabel1\n$runlabel2\n$startdate\n$enddate\n";
  print TMP "$CWFDUMP[0]\n$CWFDUMP[1]\n" if $CWF;
  close(TMP);

  foreach my $exclu ( @exclus) {
    print "starting $PROGRAM with\n".
    "EXCLU $exclu\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2\n".
    "startdate $startdate\nenddate $enddate\n";
    print "CWFDUMP1 $CWFDUMP[0]\nCWFDUMP2 $CWFDUMP[1]\n" if $CWF;

    if ($DRY_RUN) {
      print "DRY_RUN: not running '| mpirun ./$LPROG'\n";
      system("ls -lht --time-style=long-iso *")
    } else {
      open (PROG, "| mpiexec ./$LPROG") || die "Unable to execute $LPROG. Exiting.\\n" ;
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
Emissions: $emisdir
Emislist: @emislist
Meteo: $MetDir
Version: $testv
Chemical scheme: $Chem
@packages
Processors $NDX $NDY
SR?  $SR
CWF? $CWF
BC? $cwfbc
iyr_trend: $iyr_trend
------------------------------
femis: femis.$scenario
------------------------------
EOT
  close RUNLOG;

# BENCHMARK summary info
  if (%BENCHMARK and -s "$runlabel1\_fullrun.nc") {
    system("cdo infov $runlabel1\_fullrun.nc > $runlabel1\_fullrun.infov");
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
  unlink ( @list_of_files );

#tar sites and sondes. Use sondes to check as these are produced les frequently.
  my $last_sondes = sprintf  "sondes.%02d%02d", $mm2, $yy;
  print "LOOKING FOR LAST SITES $last_sondes\n";
  if ( -r $last_sondes ) {
    print "FOUND LAST sondes $last_sondes\n";
    system("tar cvzf $runlabel1.sites.tgz  sites.*");
    system("tar cvzf $runlabel1.sondes.tgz sondes.*");
  }

  if ($CWF) {
    my $old="EMEP_OUT.nc";
    my $new="$CWFDUMPDIR/$scenario\_dump.nc";    # today's dump
    system("mkdir -p $CWFDUMPDIR/; mv $old $new") if (-e "$old");
    if ($SR) {
      ($old=$new)=~s/$CWFBASE/$CWFDATE[0]/g;      # yesterday's dump
      system("rm $old") if (-e $old);
      $old="modelrun.finished";
      foreach my $task ('PBS_ARRAY_INDEX', 'PBS_ARRAYID', 'TASK_ID') {
        $new="runsr_$ENV{$task}.finished" if $ENV{$task};
        system("mkdir -p ../CWF_$CWFBASE/;echo $scenario >> ../CWF_$CWFBASE/$new")
          if (-e $old) && ($ENV{$task});
      }
    }
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

##############################################################
##### Stuff for Source receptor matrisses             ########
##############################################################
package EMEP::Sr;

my (%country_nums, @eu15, @euNew04, @eu25, @euNew06, @eu27, @sea, @noneu, @emep, @eecca, @eccomb);
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
 RFE => 74, KZE => 75, UZ => 76, TM  => 77, UZE => 78,
 TME => 79, CAS => 80, TJ => 81, ARL => 82, ARE => 83,
 ASM => 84, ASE => 85, AOE => 86,
 RFX => 87, ASX => 88, PAX => 89, AOX => 90,
 NAX => 91,
 KZT => 92, RUE => 93, UZT => 94, TMT => 95, AST => 96,
 BA2 => 302, BA3 => 303, BA4 => 304, BA5 => 305, BA6 => 306, # Baltic sep.
 BA7 => 307, BA8 => 308, BA9 => 309,
 NS2 => 312, NS3 => 313, NS4 => 314, NS5 => 315, NS6 => 316, # N. Sea sep.
 NS7 => 317, NS8 => 318, NS9 => 319,
 AT2 => 322, AT3 => 323, AT4 => 324, AT5 => 325, AT6 => 326, # Atlant. sep.
 AT7 => 327, AT8 => 328, AT9 => 329,
 ME2 => 332, ME3 => 333, ME4 => 334, ME5 => 335, ME6 => 336, # Medit. sep.
 ME7 => 337, ME8 => 338, ME9 => 339,
 BL2 => 342, BL3 => 343, BL4 => 344, BL5 => 345, BL6 => 346, # Bl. Sea sep.
 BL7 => 347, BL8 => 348, BL9 => 349, ALL => 0
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
  $base = "CWF_$CWFBASE" if $CWF ;
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
  if( $pollut eq "P" ) { $pm25 = $pmco = $redn  };
  if( $pollut eq "NP" ) { $nox = $pm25 = $pmco = $redn  };
  if( $pollut eq "SNP" ) { $sox = $nox = $pm25 = $pmco =  $redn  };
  if( $pollut eq "AN" ) { $nh3 = $nox =  $redn  };
  if( $pollut eq "SNAV" ) { $sox = $nox = $nh3 = $voc =  $redn  };
 #if( $pollut eq BASE ) then no change!

  my $femisdat = slurp("$DataDir/femis.dat");

  die "ERROR!! No country Num for $cc!\n" unless defined(my $ccnum = $country_nums{$cc});

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
  if ( $cc eq "KZT" ) {  # Add splitted countries
    foreach my $cx (53, 75) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "TMT" ) {  # Add splitted countries
    foreach my $cx (77, 79) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "UZT" ) {  # Add splitted countries
    foreach my $cx (76, 78) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "RU" ) { # Add splitted and external RU
    foreach my $cx (36..38, 42) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "RUE" ) { # Add RU and RFE
    foreach my $cx (36, 37, 38, 42, 71, 74) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "ASI" ) { # Add splitted ASI
    foreach my $cx (76, 77, 80, 82, 84) {
      $femisdat .= "$cx $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
    }
  }
  if ( $cc eq "AST" ) { # Add ASI and ASE
    foreach my $cx (80, 82, 83, 84, 85) {
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
