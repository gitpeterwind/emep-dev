#!/usr/bin/perl

#Common script for Njord, Stallo.
#Choose $VILJE=1 or $STALLO=1
#___________________________________________________________________
# queue commands for PBS

#Queue system commands start with #PBS (these are not comments!)
# Vilje: (take out one # and put one # before the Stallo). 
#   select= number of nodes, ncpus=number of threads per node to reserve, 
#   mpiprocs=number of MPI threads per node. For 64 processors:
##PBS -l select=4:ncpus=32:mpiprocs=32 -v MPI_MSGS_MAX=2097152
# Stallo:
#   Some nodes on Stallo have 16 some others 20 cpus
#   use ib for infiniband (fast interconnect).
#   lnodes= number of nodes, ppn=processor per node (max16 or 20 on stallo)
#   lpmeme=memory to reserve per processor (max 16GB per node)
#PBS -lnodes=64 -lpmem=1000MB
# Wall time limit of run
#PBS -lwalltime=00:20:00
# Make results readable for others:
#PBS -W umask=0022
# Account for billing
#PBS -A nn2890k
# Multiple tasks for paralel SR runs (one task per country)
##PBS -t 1-56
#___________________________________________________________________


#___________________________________________________________________
#Abel queue commands ?

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
foreach my $key  qw ( PBS_O_HOST HOSTNAME PBS_SERVER MACHINE ){
 next unless defined $ENV{$key};
 $STALLO = 1 if $ENV{$key} =~/stallo/;
 $VILJE  = 1 if $ENV{$key} =~/vilje/;
 #my $val = $ENV{$key};
 #print "KEY $key Val $val  STALLO?? $STALLO VILJE $VILJE\n";
}
print "HOST DETECT: V$VILJE S$STALLO \n";

# -j4 parallel make with 4 threads
my @MAKE = ("gmake", "-j4", "MACHINE=snow");
   @MAKE = ( "make", "-j4", "MACHINE=vilje")  if $VILJE==1 ;
   @MAKE = ( "make", "-j4", "MACHINE=stallo") if $STALLO==1 ;
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
my ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("2708"    ,"EmChem09soa","EMEPSTD","EMEPSTD","EECCA","EMEP");
#  ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("test"    ,"EmChem09"   ,"EMEPSTD","EMEPSTD","EECCA",0);
#  ($testv,$Chem,$exp_name,$outputs,$GRID,$MAKEMODE) = ("testcri2","CRI_v2_R5"  ,"CRITEST","EMEPSTD","EECCA",0);
#eg ($testv,$Chem,$exp_name,$GRID,$MAKEMODE) = ("tests","EmChem09","TESTS","RCA","EmChem09");

my %BENCHMARK;
# OpenSource 2008
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2005,emis=>"Modrun07/OpenSourceEmis"           ,chem=>"EmChem03");
# Dave's preference for EMEP:
#  %BENCHMARK = (grid=>"EMEP"  ,year=>2006,emis=>"Modrun10/EMEP_trend_2000-2008/2006",chem=>"EmChem09");
# EECCA Default:
## %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun11/EMEP_trend_2000-2009/2008",chem=>"EmChem09soa",make=>"EMEP");
# Status Runs:
#  %BENCHMARK = (grid=>"EECCA" ,year=>2007,emis=>"Modrun09/2009-Trend2007-CEIP") ;
#  %BENCHMARK = (grid=>"EECCA" ,year=>2008,emis=>"Modrun10/2010-Trend2008_CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2009,emis=>"Modrun11/2011-Trend2009-CEIP");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2010,emis=>"Modrun12/2012-Trend2010-CEIP",chem=>"EmChem09soa",make=>"EMEP2010");
#  %BENCHMARK = (grid=>"EECCA" ,year=>2011,emis=>"Modrun13/2013-Trend2011-CEIP",chem=>"EmChem09soa",make=>"EMEP2011");
# Alternative domains:
#  %BENCHMARK = (grid=>"TNO28" ,year=>2008,emis=>"emis_TNO28"         );
#  %BENCHMARK = (grid=>"MACC02",year=>2008,emis=>"2008_emis_EMEP_MACC") ;
if (%BENCHMARK) {
  $BENCHMARK{'archive'} = 1;                        # save summary info in $DataDir
  $BENCHMARK{'debug'} = $BENCHMARK{'archive'};      # chech if all debug flags are .false.
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

my $CWF=0;     # Set to N for 'N'-day forecast mode (0 otherwise)
   $CWF=0 if %BENCHMARK;
my ($CWFBASE, $CWFDAYS, $CWFMETV, @CWFDATE, @CWFDUMP, $eCWF, $aCWF) if $CWF;
if ($CWF) {
  $CWFBASE=$ENV{"DATE"}?$ENV{"DATE"}:"today"; # Forecast base date     (default today)
  $CWFDAYS=$ENV{"NDAY"}?$ENV{"NDAY"}:$CWF;    # Forecast lenght indays (default $CWF)
  $CWFMETV=$ENV{"UTC"}?$ENV{"UTC"}:"12";      # Met.UTC version        (default 12UTC)
  $CWFBASE=shift if @ARGV;              # Forecast base date, lenght
  $CWFDAYS=shift if @ARGV;              #  & MetUTC version can be passed
  $CWFMETV=shift if @ARGV;              #  as argument to script
  $CWFBASE="tomorrow" if($CWFBASE eq "today")and($CWFMETV =~ /12/);  # default date for 12UTC version 
  $CWFBASE=date2str($CWFBASE,"%Y%m%d");
# $CWFMETV:
# Forecast/Analysis ($CWFDAYS<=10): meteo${CWFBASE}_{00,01,..}d.nc
#   AN00|AN Analysis w/ 00UTC met
#   FC12|12 Forecast w/ 12UTC met
# Hindcast ($CWFDAYS>10): meteo{DAY1,DAY2,..}_??d.nc
#   *00|24|48|72 run w/ 00UTC met 00d|01d|02d|03d
#   *12|36|60|84 run w/ 12UTC met 01d|02d|03d|04d
  $eCWF=0;                              # Emergency forecast
  $aCWF=($CWFMETV =~ /AN/ );            # Analysis
  $CWF=($eCWF?"eemep-":"CWF_").($CWFMETV?"$CWFMETV-$CWFBASE":"$CWFBASE");
  $CWFMETV =~s/[^\d.]//g;                           # extract number part
  $CWFDATE[0]=date2str($CWFBASE." 1 day ago"  ,"%Y%m%d");     # yesterday
  $CWFDATE[1]=$CWFBASE;                                       # start date
  $CWFDATE[2]=date2str($CWFDATE[0]." $CWFDAYS day","%Y%m%d"); # end date
##$CWFDUMP[0]=date2str($CWFBASE."                 ,"%Y-1-1"); # dump/nest every day at 00
  $CWFDUMP[0]=date2str($CWFBASE." 1 day"          ,"%Y%m%d"); # 1st dump/nest
  $CWFDUMP[1]=date2str($CWFBASE." 2 day"          ,"%Y%m%d"); # 2nd dump/nest
  $MAKEMODE=($eCWF)?"eEMEP":"MACC";    # Standard Forecast model setup
  $MAKEMODE .="-3DVar" if($aCWF);
  $exp_name = ($eCWF)?"EMERGENCY":($aCWF?"ANALYSIS":"FORECAST");
  $testv.= ($eCWF)?".eCWF":".CWF";
##$MAKEMODE=$eCWF?"eEMEP2010":"MACC-EVA2010";    # 2010 special
# $MAKEMODE=$eCWF?"eEMEP2011":"MACC-EVA2011";    # 2011 special
# $exp_name.=_MACCEVA;
##$MAKEMODE=$eCWF?"eEMEP2013":"MACC";            # 2013 special
  $GRID = ($eCWF)?"GLOBAL":"MACC14";
}
 $MAKEMODE="SR-$MAKEMODE" if($MAKEMODE and $SR);

# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables
#      are explained below, and derived variables set later.-

my $year = "2006";
   $year = substr($CWFBASE,0,4) if $CWF;
   $year = $BENCHMARK{"year"} if %BENCHMARK;
( my $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right

# iyr_trend:
# :can be set to meteorology year or arbitrary year, say 2050

my $iyr_trend = $year;
#$iyr_trend = "1990" ;  #  RCA ECLAIRE

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

my $MONTHLY_EMIS = ($GRID eq "GLOBAL"); #Switch off if only annual used
   $MONTHLY_EMIS = 0 if ($exp_name =~ /ECLAIRE/);

my ($HOMEROOT, $WORKROOT, $MetDir);
our $DataDir;
if ($STALLO) {
  $HOMEROOT = "/home";
  $WORKROOT = "/global/work";
  $DataDir  = "$WORKROOT/$PETER/Data";
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year";
  $MetDir   = "$DataDir/$GRID/metdata_H20/$year" unless -d $MetDir;
  $MetDir   = "$DataDir/$GRID/metdata/$year"     unless -d $MetDir;
  $MetDir   = "$WORKROOT/$PETER/ClimData/$year"  if ($GRID eq "RCA" );
  $MetDir   = "$DataDir/$GRID/metdata_CWF/$year" if $CWF;
} else { #Ve or Vilje
  $HOMEROOT = "/home/metno";
  $WORKROOT = "$HOMEROOT/$USER/work";
  $DataDir  = "$HOMEROOT/mifapw/work/Data";
  $MetDir   = "$DataDir/$GRID/metdata_EC/$year";
  $MetDir   = "$DataDir/$GRID/metdata_CWF/$year" if $CWF;
  $MetDir   = "/prod/forecast/work/emep/ec/prepmet" if $eCWF and ($USER eq $FORCAST);
}
#$MetDir   = "/hard/code/path/to/$GRID/metdata/$year/if/necessary";
die "Missing MetDir='$MetDir'\n" unless -d $MetDir;

# DataDir    = Main general Data directory
my $DATA_LOCAL = "$DataDir/$GRID";   # Grid specific data , EMEP, EECCA, GLOBAL
# Pollen data
my $PollenDir = "$HOMEROOT/$BIRTHE/Unify/MyData";
   $PollenDir = 0 unless -d $PollenDir;
# Eruption (eEMEP)
my $EmergencyData = "$HOMEROOT/$ALVARO/Unify/MyData";
   $EmergencyData = 0 unless -d $EmergencyData;
# Project-specific data:
my $ProjDataDir = "";          # Change for specific project data
#  $ProjDataDir = "$WORKROOT/mifads/Data/inputs_projects/eclaire_Jun2013"; #e.g.


# Boundary conditions: set source direcories here:
# BCs can come from Logan, Fortuin, UiO (CTM2) or EMEP model runs:

my $VBS   = 0;

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

my $WORKDIR = "$WORKROOT/$USER/$testv.$year";  # working and result directory
   $WORKDIR = "$WORKROOT/$testv.$year" if($WORKROOT =~ /$USER/);
   $WORKDIR =~ s|$testv.$year|Benchmark/$GRID.$year|g if (%BENCHMARK);
   $WORKDIR = "/prod/forecast/run/eemep" if $eCWF and ($USER eq $FORCAST);
my $MyDataDir   = "$HOMEROOT/$USER/Unify/MyData";           # for each user's private input
my $SoilDir     = "$DATA_LOCAL/dust_input";               # Saharan BIC
   $SoilDir = 0 unless -d "$SoilDir/BC_DUST/2000";

# TEST! Road dust NOTE! The road dust code may not be working properly yet! Not tested enough!
my $RoadDir = "$HOMEROOT/$ROBERT/Unify/MyData/TNO_traffic" ;
   $RoadDir = "$HOMEROOT/$AGNES/MyData/TNO_traffic" if $VILJE;
   $RoadDir = 0 unless -d $RoadDir;
   $RoadDir = 0 if $CWF;

# Forecast: nest/dump dir, BCs pattern
my ($CWFIC, $CWFBC, $CWFPL) if $CWF;
my ($cwfic, $cwfbc, $cwfpl) = ("No IC file","No BC file","No Pollen file");
if ($CWF) {
 ($CWFIC  = "${CWF}_dump.nc" ) =~ s|$CWFBASE|%Y%m%d|g;
 ($CWFIC  = "$WORKDIR/$CWFIC") =~ s|$testv.$year|$testv.dump|;
  $CWFIC  =~ s|run/eemep|work/emep/restart| if $eCWF and ($USER eq $FORCAST);
  $CWFBC  = "$DataDir/$GRID/Boundary_conditions/";            # IFS-MOZ
  $CWFBC .= "%Y_IFS-MOZART_FC/cwf-mozifs_h%Y%m%d00_raqbc.nc";# :Forecast
# $CWFBC .= "%Y_IFS-MOZART_AN/h%Y%m%d00_raqbc.nc";           # :ReAnalysus
# $CWFBC .= "%Y_EVA/EVA_%Y%m%d_EU_AQ.nc";                    # :EVA-2010/2011
 ($CWFPL  = $CWFIC) =~ s|_dump|_pollen|;
}

#ds check: and change
chdir "$ProgDir";
#die "Dir wrong!!!!! $testv label does not match in ENV$ENV{PWD}\n"
#  unless $ENV{PWD} =~ /Unimod.$testv.$year/;
print "TESTING ENV:", $ENV{PWD}, "\n";


# Default emissplits used here. if $Specials is set will look
my $SplitDir = "$DataDir/SPLITS_JAN2010/BASE_NAEI2000_GH2009.$Chem" ;
   $SplitDir = "$ChemDir/EMISSPLIT";
#RB:had "~mifarb/Unify/MyData/D_EGU/SPLITS_NOV2009v2/BASE_NAEI2000_GH2009.$Chem" ;

my $version     = "Unimod" ;
   $version     = "ZD_3DVar/Unimod_3DVar" if($aCWF);
my $PROGRAM     = "$ProgDir/$version";        # programme
my $subv        = "$testv" ;                  # sub-version (to track changes)

# New system. For "normal runs", we just put one run-name into the array @runs.
# For SR runs we can add many scenarios - dealt with later.
# The effect is to choose the approproate femis file

my $scenario = "Base";     # Reset later if SR
my @runs     = ( $scenario );

#Possible emission scenarios for HIRHAM run
#GEA scenarios: HIGH_CLE, HIGH_FROZEN, LOW_CLE, LOW_SLE
#Historic emissions from Lamargue et al. : historic_emis
my ($emisscen,$emisyear) = ("historic_emis","$year");
    $scenario="${emisscen}_emis${emisyear}_met${year}" if $GRID eq "HIRHAM";

#EMISSIONS: default settings
my ($EMIS_INP,$emisdir,$pm_emisdir)=("$DATA_LOCAL","none","none");

given($GRID){
  when("EECCA"){given($year){
    when([2000..2009]){$emisdir="$EMIS_INP/Modrun12/EMEP_trend_2000-2009/$year";}
    when([2010])      {$emisdir="$EMIS_INP/Modrun12/2012-Trend$year-CEIP";}
    when([2011])      {$emisdir="$EMIS_INP/Modrun13/2013-Trend$year-CEIP";}
  }}
  when("EMEP"){given($year){
    when([1990..1999]){$emisdir="/global/work/$AGNES/Emission_Trends/$year";}
    when([2000..2008]){$emisdir="$EMIS_INP/Modrun11/EMEP_trend_2000-2009/$year";}
    when([2009])      {$emisdir="$EMIS_INP/Modrun11/2011-Trend$year-CEIP";}
  }}
  when(/TNO/){
    $emisdir="$WORKROOT/$AGNES/emis_SRbase/INERIS_direct/$GRID"      if($STALLO);
    $emisdir="$HOMEROOT/$AGNES/emission/SD_emis/INERIS_direct/$GRID" if($VILJE);
  }
  when("GLOBAL"){
    $emisdir=($eCWF)?"$EMIS_INP/Emissions_June2012":"$EMIS_INP/MonthlyEmis";
    $emisdir="$EMIS_INP/ECLAIRE_1deg_Feb2013/$iyr_trend" if($exp_name=~/ECLAIRE/);
  }
  when("RCA")    {$emisdir="$ProjDataDir/Interpolations";} #EnsClim
  when("HIRHAM") {$emisdir="$EMIS_INP/emissions/$emisscen/$emisyear";}
  when("GEMS025"){$emisdir="$EMIS_INP/Emissions/2008-Trend2006-V9-Extended_PM_corrected-V3";}
# when("MACC02") {$emisdir="$EMIS_INP/Emissions/2008_emis_EMEP_from_PS50";}
# when("MACC02") {$emisdir="$EMIS_INP/Emissions/2008_emis_EMEP_MACC";}
  when("MACC02") {$emisdir="$EMIS_INP/Emissions/2007_emis_MACC";}
  when("MACC14") {$emisdir="$EMIS_INP/Emissions/2009_emis_MACCII";}
}
#TMP and should be improved because it gives errors for other domains!
#.. For using emissions of EC/OC instead of PMx
my $RFEmisDir  = "/global/work/$SVETLANA/Data_RF";  # Split-Fraction files for EC/OC
my $TNOemisDir = "/global/work/$SVETLANA/Emis_TNO"; # TNO EC/OC emissions
#$emisdir = $TNOemisDir if $EUCAARI;

given($GRID){
  when("GEMS025"){ $pm_emisdir = $emisdir; }
  when("MACC02",){ $pm_emisdir = $emisdir; }
  when("MACC14",){ $pm_emisdir = $emisdir; }
  when(/TNO/)    { $pm_emisdir = $emisdir; }
  when("GLOBAL") { $pm_emisdir = $emisdir; }
  default {
    $pm_emisdir = $emisdir;
    $pm_emisdir = "$EMIS_INP/2006-Trend2000-V7"  if $year < 2000;
    $pm_emisdir = "/home/$ROBERT/Unify/MyData/D_EGU/${GRID}_GRID" if $VBS;
  }
}

#EMISSIONS: BENCHMARK settings
if (%BENCHMARK){
  $emisdir    = "$EMIS_INP/$BENCHMARK{'emis'}";
  $pm_emisdir = $emisdir;
}
my $CDF_EMIS=0;#put one if TNO7.nc emissions are used
die "Missing emisdir='$emisdir' for GRID='$GRID'\n"       unless (-d $emisdir or $CDF_EMIS);
die "Missing pm_emisdir='$pm_emisdir' for GRID='$GRID'\n" unless (-d $pm_emisdir or $CDF_EMIS);

#FEB 2013 TEST of netcdf emissions
my $SNAP_CDF = "/global/work/mifapw/temp";  # Use for CdfFractions
   $SNAP_CDF = "/global/work/mifads/cdf_emis";  # Use for CdfSnap
   $SNAP_CDF = 0 unless $STALLO;

#Dave, reset to Emission_Trends for Chem project, Oct 18th
my $TREND_RUNS = 0;
if ($STALLO && $TREND_RUNS ) {
  $EMIS_INP = "/global/work/$AGNES/Emission_Trends";
  die "Year not in trend run series!! " unless -f $EMIS_INP/$year;
  $emisdir = "$EMIS_INP/$year";
  $pm_emisdir = $emisdir;
}

my $RESET        = 0 ;  # usually 0 (false) is ok, but set to 1 for full restart
my $COMPILE_ONLY = 0 ;  # usually 0 (false) is ok, but set to 1 for compile-only
my $DRY_RUN      = 0 ;  # Test script without running model (but compiling)

if (%BENCHMARK and $BENCHMARK{'debug'}){
  die "No debug flags for benchmarks!"
  if (system ("grep -Hnie 'logical.*DEBUG.*=\ *.TRUE.' $ProgDir/*.f90") == 0) or
     (system ("grep -Hnie 'DEBUG.*=\ *.TRUE.' $ProgDir/ModelConstants_ml.f90") == 0);
}

if ( $ENV{PBS_NODEFILE} ) {
   $_ =  `wc -l $ENV{PBS_NODEFILE}`;
   my $RUN_NPROC;
   ( $RUN_NPROC, undef ) = split;
   print "Qsub required: $RUN_NPROC processors\n";

} else {
   print "skip nodefile check on interactive runs\n";
}


my @month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

#Only 360 days in HIRHAM metdata. We ignore leaps
@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,24) if $GRID eq "HIRHAM";

my $mm1 ="06";      # first month, use 2-digits!
my $mm2 ="06";      # last month, use 2-digits!
my $dd1 =  1;       # Start day, usually 1
my $dd2 =  0;       # End day (can be too large; will be limited to max number of days in the month)
                    # put dd2=0 for 3 hours run/test.
# Allways runn full year on benchmark mode
($mm1,$mm2,$dd1,$dd2)=("01","12",1,31) if (%BENCHMARK);

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


if ( $RESET ) { ########## Recompile everything!

  # For now, we simply recompile everything!
  system(@MAKE, "clean");
  if ($MAKEMODE) {
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

if ($SR or defined($ENV{'PBS_ARRAY_INDEX'}) or 
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
  my $runlabel2 = "$testv\_$Chem\_$scenario\_$year\_Trend$iyr_trend";   # NO SPACES! LONG (written into CDF files)

  my $RESDIR = "$WORKDIR/$scenario";
     $RESDIR = "$WORKDIR/$scenario.$iyr_trend" if ($GRID eq "RCA");
  mkdir_p($RESDIR);

  chdir $RESDIR;   ############ ------ Change to RESDIR
  print "Working in directory: $RESDIR\n";

  # Meteorology name in felt/cdf format
  $METformat=($METformat eq "felt")?"$MetDir/fhh.YYYYMMDD":      # felt
                                    "$MetDir/meteoYYYYMMDD.nc";  # cdf
  if($CWF) { # Forecast Meteorology in NetCDF
    $METformat="./meteoYYYYMMDD.nc";            # link file to work path
    my $metday = 0;
    if($CWFMETV) { # UTC version and DAY offset for MET.UTC version
      $metday = int($CWFMETV/24+0.99);             # DAY offset
      $MetDir.= sprintf("_%02dUTC",($CWFMETV%24)) # 00/12 UTC versions
                unless ($USER eq $FORCAST);
      die "Missing MetDir='$MetDir'\n" unless -d $MetDir;
    }
    $MetDir=~s:$year:%Y:g;                      # Genereal case for Jan 1st
    for (my $n = 0; $n < $CWFDAYS; $n++) {
      my $metfile="$MetDir/meteo%Y%m%d_%%02d.nc";
      if($CWFDAYS<=10){  # forecast with 00/12 UTC FC met.
        $metfile=sprintf date2str($CWFBASE." $metday day ago",$metfile),$n+$metday;
      }else{             # hindcast using day 0,..,3 FC met.
        $metfile=sprintf date2str($CWFBASE." $metday day ago $n day",$metfile),$metday;
      }
      # Chech if meteo is in place
      die "Meteo file $metfile for $CWFBASE not available (yet). Try later...\n"
        unless (-e $metfile);
      $CWFDATE[2]=date2str($CWFBASE." $n day","%Y%m%d");
      mylink("CWF Met:",$metfile,date2str($CWFDATE[2],$METformat)) ;
      # IFS-MOZART BC file
      $cwfbc=date2str($CWFDATE[2],$CWFBC);
      $cwfbc=date2str($CWFDATE[0],$CWFBC) unless (-e $cwfbc);
      mylink("CWF BC :",$cwfbc,"EMEP_IN_BC_${CWFDATE[2]}.nc") if (-e $cwfbc);
    }
# Forecast nest/dump files
    $cwfic=date2str($CWFDATE[0],$CWFIC);  # yesterday's dump
    if($aCWF){                              # CWF_AN run: if dump not found
      $cwfic=~s/AN-/FC-/ unless(-e $cwfic); # use CWF_FC dump
      $cwfic=~s/FC-//    unless(-e $cwfic); # use CWF dump
      die "$CWF restart file for $CWFBASE not available (yet):\n\t$CWFIC\n"
          ."Try later...\n" unless(-e $cwfic);
    }else{                                  # CWF_FC run: always from Analysis
      $cwfic=~s/FC-/AN-/;                   # use CWF_AN dump
    }
    my $metfile="$MetDir/meteo${CWFDATE[0]}_00.nc"; # yesterday's met
    if (-e $metfile and ! -e $cwfic) {
      print "No dumpfile present:\n\t$cwfbc.\n\tTake extra spin-up day.\n";
      # Update simulation: start-date ($CWFDATE[1]), "yesterday" ($CWFDATE[0])
      $CWFDATE[1]=$CWFDATE[0];
      $CWFDATE[0]=date2str($CWFDATE[0]." 1 day ago","%Y%m%d");
      mylink("CWF Met:",$metfile,date2str($CWFDATE[1],$METformat));
      # IFS-MOZART BC file
      $cwfbc=date2str($CWFDATE[1],$CWFBC);
      $cwfbc=date2str($CWFDATE[0],$CWFBC) unless (-e $cwfbc);
      mylink("CWF BC :",$cwfbc,"EMEP_IN_BC_${CWFDATE[1]}.nc") if (-e $cwfbc);
      # see if we can link to a dump file ...
      $cwfic=date2str($CWFDATE[0],$CWFIC);  # yesterday's dump
     #$cwfic=~s/AN-//        if($aCWF);     # use Forecast dump on AN runs
      $cwfic=~s/FC-/AN-/ unless($aCWF);     # use Analyis dump on FC runs
    }
    mylink("CWF IC :",$cwfic,"EMEP_IN_IC.nc") if (-e $cwfic);
# Update start and end months (used for linking some climatological files)
    $mm1=substr($CWFDATE[1],4,2);  # start date
    $mm2=substr($CWFDATE[2],4,2);  # end date
  }

#Experimental: Grid definitions in a separate file 
  my $old = "$DATA_LOCAL/Grid_Def.nc";
  my $new = "Grid_Def.nc";
  mylink( "Linking:", $old, $new);
#  my $old = "$ProgDir/Vertical_levels20.txt";
#  my $new = "Vertical_levels.txt";
#  mylink( "Linking:", $old, $new);

#=================== INPUT FILES =========================================
# ToDo Change noxsplit.default to defaults, as with voc (also in Unimod)

  my %ifile = ();   # List of input data-files
  my %inml  = ();   # List of input data-files in nml format

# First, emission files are labelled e.g. gridSOx, which we assign to
# emislist.sox to ensure compatability with the names (sox,...) used
# in the model. It doesn't matter if we have extra mapping here.
# GenIn.reactions and the associated emislist decides what gets used.
# e.g. lines such as:
#  emisfiles:sox,nox,co,voc,nh3
#  emisfiles:pm25
# etc.

  my $timeseries  = "$DataDir/inputs_emepdefaults_Jun2012";
  my $Tbase = 18 ;  # Base-temperature for Degree-day (HDD) files
  ($timeseries,$Tbase)=("$DataDir/inputs_eurodelta_Jun2012",20) if($INERIS_FACS);

  #$ifile{"$timeseries/HourlyFacs.EMEP2003"} = "HOURLY-FACS";
  #$ifile{"$timeseries/HourlyFacs.TNO2005"} = "HOURLY-FACS";
# INERIS provided the most complete hourly file, we use as default
  $ifile{"$timeseries/HourlyFacs.INERIS"} = "HOURLY-FACS";
  $ifile{"$timeseries/EmisHeights.txt"} = "EmisHeights.txt";
#  $ifile{"$timeseries/EmisHeights_P.txt"} = "EmisHeights.txt"; #not much tested yet

  my %gridmap=(
    "co"  =>"CO"  ,"nh3" =>"NH3" ,"voc"=>"NMVOC","sox" =>"SOx","nox"=>"NOx" ,
    "pm10"=>"PM10","pm25"=>"PM25","pmco"=>"PMco",
 # VBS specials
    "pocfwd"=>"POCfWD","pocffl"=>"POCfFL","poccfl"=>"POCcFL",
    "ecfwd" =>"ECfWD" ,"eccwd" =>"ECcWD" ,"ecffl" =>"ECfFL" ,
    "eccfl" =>"ECcFL" ,"forfbc"=>"FORFBC","forfoc"=>"FORFOC",
 #  Sometimes used also:
    "ecfi"  =>"ECfine","ecco"  =>"ECcoar","ocfi"  =>"OCfine");
                         # sometimes was  "ocfi" =>"POCfine");

  foreach my $poll (@emislist) {
    my $dir = $emisdir;
    $dir = $pm_emisdir if $poll =~ /pm/;   # FIX needed prior to 2000
 # VBS specials #rb Wood burning, Fossil fuel and Forest fire PM from TNO files
    $dir = $pm_emisdir if $poll =~ /wd/;   #
    $dir = $pm_emisdir if $poll =~ /fl/;   #
    $dir = $pm_emisdir if $poll =~ /forf/;   #
    print "TESTING PM $poll $dir\n";

    if ($GRID eq "MACC14") { # For most cases only Emis_TNO7.nc is available
      $ifile{"$emisdir/EmisOutFrac.$poll"} = "emislist.$poll" if(-e "$emisdir/EmisOutFrac.$poll");
    }elsif( $GRID eq "RCA"){
      #EnsClim RCA #$ifile{"$dir/grid$gridmap{$poll}"} = "emislist.$poll";
      #$ifile{"$emisdir/EmisOut_2005.$poll"} = "emislist.$poll";
      #$ifile{"$emisdir/EmisOut_$iyr_trend.$poll"} = "emislist.$poll";
      # ECLAIRE - uses 2005 as base for historical
      $ifile{"$emisdir/Emis_RCA_Grid_2005/EmisOut_2005.$poll"} = "emislist.$poll";
    }elsif(($NH3EMIS_VAR)&&($poll eq "nh3")){
      $dir = "$HOMEROOT/$AGNES/emis_NMR";
      $ifile{"$dir/gridNH3_NMR_$year"} = "emislist.$poll";
    }elsif($CDF_EMIS){
	#no linking
    }else{
      $ifile{"$dir/grid$gridmap{$poll}"} = "emislist.$poll";
    }
    $dir=(-e "$emisdir/Emis_TNO7.nc")?$emisdir:$DataDir;
    $ifile{"$dir/Emis_TNO7.nc"} = "EmisFracs_TNO7.nc";

    if($SNAP_CDF) { # in testing:
      print "SNAP CDF TESTS $poll\n";
      #$ifile{"$SNAP_CDF/Emis_$gridmap{$poll}.nc"}
      #2005:
      $ifile{"$SNAP_CDF/MACC2_Mar2013/2005/Emis_$gridmap{$poll}.nc"}
                        = "GriddedSnapEmis_$poll.nc" if $SNAP_CDF ;

      #2005 $ifile{"$SNAP_CDF/CdfGlobal/$iyr_trend/Emis_$gridmap{$poll}.nc"}  #2474erca had RCAmap?
      $ifile{"$SNAP_CDF/CdfGlobal/2005/Emis_$gridmap{$poll}.nc"}  #2474erca had RCAmap?
         = "GlobalSnapEmis_$poll.nc" if $SNAP_CDF ;

      my $ship = "$SNAP_CDF/CdfGlobal/IPCC_v1_20_04_2009_emep/EmisIPCC_$gridmap{$poll}_ships_$iyr_trend.nc";
      $ship = "$SNAP_CDF/CdfGlobal/IPCC_v1_20_04_2009_emep/EmisIPCC_$gridmap{$poll}_ships_2000.nc"; #TESTING with 2000
      $ifile{$ship} = "GlobalShipEmis_$poll.nc" if -e $ship ;
      #die "SNAP TEST WILL NOT WORK FOR $ship THIS YEAR$iyr_trend  \n" unless -f $ship; #Only 1990, 1990, steps of 10 so far
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

    if ( $Chem eq "EmChem09" ) { # e.g. when PM25 is not split, e.g. RCA, make EMCHEM09
      $ifile{"$SplitDir/emissplit.specials.$poll"} = "emissplit.specials.$poll"
      if( -e "$SplitDir/emissplit.specials.$poll" );
    } elsif ( $Chem eq "CRI_v2_R5" ) { # e.g. TSAP
       print "NO SPECIALS in EMISSPLIT for $Chem DIR was $SplitDir\n";
    } elsif ( -e "$timeseries/emissplit.$Specials.$poll.$iyr_trend" ) { # e.g. TSAP
      $ifile{"$timeseries/emissplit.$Specials.$poll.$iyr_trend"} =
             "emissplit.specials.$poll"
    } elsif ( -e "$timeseries/emissplit.specials.$poll.2010" ) { # when no other year availanle
      $ifile{"$timeseries/emissplit.specials.$poll.2010"} =
             "emissplit.specials.$poll"
    } elsif ( -e "$SplitDir/emissplit.specials.$poll" ) {
      $ifile{"$SplitDir/emissplit.specials.$poll"} =
             "emissplit.specials.$poll";
    }
  }

  foreach my $mmm ( $mm1 .. $mm2, $mm1, $mm2 ) {
    my $mm = sprintf "%2.2d", $mmm;
    $ifile{"$DATA_LOCAL/natso2$mm.dat"} =  "natso2$mm.dat" unless ($GRID eq "MACC14");
    $ifile{"$DataDir/lt21-nox.dat$mm"} =  "lightning$mm.dat";
# BIC for Saharan dust
    if ( $SoilDir ) { # Not yet for EMEP domain
      foreach my $bc ( qw ( DUST_c_ext DUST_f_ext )) { #
        $ifile{"$SoilDir/BC_DUST/2000/$bc.$mm"} =  "$bc.$mm";
      }
    } # dust
    if($GRID eq "GLOBAL" && $MONTHLY_EMIS) {
      $mm = sprintf "%2.2d", $mmm;
      foreach my $t ( qw (nox voc co nh3 pm25 pmco) ) {
        $ifile{"$emisdir/grid$gridmap{$t}.$mm"} =  "grid$t.$mm";
      }
      $ifile{"$emisdir/gridSO2.$mm"} =  "gridsox.$mm";
    }
  }

# Emissions setup:
#  if ($EUCAARI) { # DS RE-CHECK shouldn't be needed
#    $ifile{"$TNOemisDir/femis.dat"} =  "femis.dat";
#    $ifile{"$DATA_LOCAL/emissions/femis.dat"} =  "femis.dat" if $GRID eq "HIRHAM" ;

  if ($ProjDataDir =~ /eclaire/ ) { # As example
    $ifile{"$ProjDataDir/femis.ecl2005to$iyr_trend"} =  "femis.dat"; 
  } else {
    $ifile{"$ChemDir/femis.defaults"} =  "femis.dat";  # created now by GenChem
  }

# my $old="$DATA_LOCAL/Boundary_and_Initial_Conditions.nc";
# my $new="Boundary_and_Initial_Conditions.nc";
# mylink( "BIC: ", $old,$new ) ;
#EUCAARI, but all?
# Skip:  $ifile{"$DATA_LOCAL/Boundary_and_Initial_Conditions.nc"} =
#                     "Boundary_and_Initial_Conditions.nc" unless ($GRID =~ /MACC/);
  $ifile{"$DataDir/Logan_P.nc"} = "Logan_P.nc";#instead of GLOBAL_O3.nc
  $ifile{"$DataDir/GLOBAL_O3.nc"} = "GLOBAL_O3.nc";
  $ifile{"$DataDir/amilt42-nox.dat"} = "ancatmil.dat";#RENAME TO AIRCARAFT?!
  $ifile{"$DataDir/GLOBAL_ForestFireEmis.nc"} =                     # GFED emissions
    "GFED_ForestFireEmis.nc";     #if ($year >= 2001 and $year <= 2007);
  $ifile{"$DataDir/ForestFire/FINN/ForestFire_Emis_$year.nc"} =     # FINN emissions
    "FINN_ForestFireEmis_$year.nc" if ($year >= 2002 and $year <= 2011);
  $ifile{"$DataDir/ForestFire/GFAS/GFAS_ForestFireEmis_$year.nc"} = # GFAS emissions
    "GFAS_ForestFireEmis_$year.nc" if ($year >= 2008 and $year <= 2012);
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

# April 2013. Now use EU emissions as proxy for future changes
#if ( $iyr_trend > 2015 )  {
#  $ifile{"$DataDir/AnnualNdep_TNO28_2020.nc"} = "annualNdep.nc";
#} else { # Note PS50x - hand-edited version
  $ifile{"$DataDir/AnnualNdep_PS50x_EECCA2005_2009.nc"} = "annualNdep.nc";
#}

# hb NH3emis
# New ammonia emissions  ---   NB no read permissions yet!!
  $ifile{"/home/$HALDIS/Unimod_NMR_NH3/Unimod.rv3_6_8/Sector_NH3Emis.txt"}="Sector_NH3Emis.txt" if($NH3EMIS_VAR);

# new inputs style (Aug 2007)  with compulsory headers:
# From rv3_14 used only for FORECAST mode
  $ifile{"$DATA_LOCAL/Inputs.Landuse"} = "Inputs.Landuse" if ($CWF and ($GRID ne "MACC14")) ;
  $ifile{"$DataDir/Landuse/landuseGLC2000_INT1.nc"} ="GLOBAL_landuse.nc";

  $ifile{"$DataDir/LanduseGLC.nc"} ="LanduseGLC.nc";
  # NB: a 1km Landuse is also available 
  $ifile{"$DataDir/Landuse/Landuse_PS_5km_LC.nc"} ="Landuse_PS_5km_LC.nc";
#  $ifile{"$DataDir/Landuse/Landuse_PS_1km_LC.nc"} ="Landuse_PS_5km_LC.nc";

  $ifile{"$DataDir/LandInputs_Mar2013/Inputs_DO3SE.csv"} = "Inputs_DO3SE.csv";
  $ifile{"$DataDir/LandInputs_Mar2013/Inputs_LandDefs.csv"} = "Inputs_LandDefs.csv";

#For dust: clay and sand fractions
  $ifile{"$DataDir/Soil_Tegen.nc"} ="Soil_Tegen.nc";

  $ifile{"$DataDir/sondesLL.dat"} = "sondes.dat";
  $ifile{"$DataDir/sitesLL.dat"} = "sites.dat";
  #DS RCA:$ifile{"$MyDataDir/sondesLLBC.dat"} = "sondes.dat";
  # Extended to get isoprene, HCHO EC, OC /(huge list!)
  #$ifile{"$MyDataDir/sitesCPM_ds.dat"} = "sites.dat";

  #LPS: point sources can  be added if needed.
  $ifile{"$MyDataDir/PointSources.txt"} = "PointSources.txt" 
   if(-e "$MyDataDir/PointSources.txt");

# DEGREE DAYS (Tbase set above, either 18 or 20):
  unless ($CWF or ($GRID eq "RCA") or ($GRID eq "GLOBAL")) {
    my $HDD = "$MetDir/HDD${Tbase}-${GRID}-$year.nc";
    print "Looking for DegreeDayFac: $HDD \n";
    system("ls -lh --time-style=long-iso $HDD");
    die "NO HDD files " unless -f $HDD;   # Can comment out if USE_DEGREEDAYS
                                          # set false in ModelConstants_ml
    $ifile{"$HDD"} = "DegreeDayFactors.nc" if -f $HDD ;
  }

  $ifile{"$DataDir/GLOBAL_LAInBVOC.nc"} = "GLOBAL_LAInBVOC.nc"; 
  #New EURO BVOC
  $ifile{"$DataDir/LandInputs_Mar2011/EMEP_EuroBVOC.nc"} = "EMEP_EuroBVOC.nc";

# Seasonal stuff  ----    Can't we improve this? e.g. every month?
  my %seasons = ("jan"=>"01","apr"=>"02","jul"=>"03","oct"=>"04");
  foreach my $s (keys %seasons) {
    $ifile{"$DataDir/a${s}t42-nox.dat"}= "ancat$seasons{$s}.dat";
    $ifile{"$DataDir/jclear.$s"}      = "jclear$seasons{$s}.dat";
    $ifile{"$DataDir/jcl1.$s"}        = "jcl1km$seasons{$s}.dat";
    $ifile{"$DataDir/jcl3.$s"}        = "jcl3km$seasons{$s}.dat";
  }

 #EnsClim RCA, and should be default:
  if ($GRID eq "RCA") {
    $ifile{"$DataDir/VolcanoesLL_2010.dat"} = "VolcanoesLL.dat";
  } else {
    $ifile{"$DataDir/VolcanoesLL.dat"} = "VolcanoesLL.dat";
  }

# Emergency senarios (eEMEP)
  if(($MAKEMODE =~ /(2010|2011)/) or ($MAKEMODE =~ /eEMEP/)){
    my $dir="$ProgDir/ZCM_Emergency";
    cp ("$dir/emergency_emission.csv","emergency_emission.csv");
    $ifile{"$dir/emergency_location.csv"} = "emergency_location.csv";
    print "$dir/emergency_location.csv\n";
    open(IN,"<$dir/emergency_location.csv");
    while(my $line = <IN>){
      unless ($line =~ /#.*/) {             # Skip comment lines
        my $vname = (split(",",$line))[0];  # Emergency tracer name
        my $efile = "$EmergencyData/${vname}_7bin.eruptions";  # Volcanic eruption
        $efile = "$EmergencyData/${vname}_2bin_${MAKEMODE}.eruptions" if($MAKEMODE =~ /(2010|2011)/);
      # $efile = "$EmergencyData/${vname}_2bin_${MAKEMODE}_SR-SOx.eruptions" if($MAKEMODE =~ /(2010|2011)/);
      # $efile = "$EmergencyData/${vname}_2bin_${MAKEMODE}_SR-PMx.eruptions" if($MAKEMODE =~ /(2010|2011)/);
        $efile =~ s|bin_SR-|bin_| if($SR); # remove the SR- from $MAKEMODE
        system("cat $efile >> emergency_emission.csv") if (-e $efile);
        $efile = "$EmergencyData/${vname}.accident";           # NPP accident
        system("cat $efile >> emergency_emission.csv") if (-e $efile);
        $efile = "$EmergencyData/${vname}.explosion";          # NUC explosion
        system("cat $efile >> emergency_emission.csv") if (-e $efile);
      }
    }
    close(IN);
  }

# For Pollen
  if($PollenDir) {
    $ifile{"$PollenDir/pollen_data.nc"} = "pollen_data.nc";
    $cwfpl=date2str($CWFDATE[0],$CWFPL) if($CWF);
    $ifile{"$cwfpl"} = "POLLEN_IN.nc"   if(-e $cwfpl);
  }

# For windblown dust
  if($SoilDir) {
    $ifile{"$SoilDir/clay_isric_percent_ext.dat"} = "clay_frac.dat";
    $ifile{"$SoilDir/sand_isric_percent_ext.dat"} = "sand_frac.dat";
  }

# TEST!!! Road dust NOTE! The road dust code is NOT thoroughly tested yet!
# NOTE ALSO THAT the Climate factors in the file below are just rough estimates based on the TNO soil water data, to be updated with something based on EMEP soil water!
  if (($GRID eq "EECCA")and($RoadDir)) {
    $ifile{"$RoadDir/RoadDust_HIGHWAYplus_emis_potential.txt"} = "HIGHWAYplus";
    $ifile{"$RoadDir/RoadDust_NonHighway_emis_potential.txt"} = "NONHIGHWAY";
    $ifile{"$RoadDir/RoughTestClimateFactorSoilWater.txt"} = "ROADDUST_CLIMATE_FAC";
  }elsif($GRID =~ /TNO/){
    $ifile{"$DATA_LOCAL/nonHIGHWAYs_RoadDust_potentials.txt"} = "NONHIGHWAY";
    $ifile{"$DATA_LOCAL/HIGHWAYplus_RoadDust_potentials.txt"} = "HIGHWAYplus";
    $ifile{"$DATA_LOCAL/ClimateFactors_SMI.txt"} = "ROADDUST_CLIMATE_FAC";
  }
# IFZ-MOZ BCs levels description (in cdo zaxisdes/eta format)
  $inml{'filename_eta'}= "$DataDir/$GRID/Boundary_conditions/mozart_eta.zaxis"
    if ( $CWF and -e $cwfbc and $cwfbc =~ m/IFS-MOZART/ );

  foreach my $f (sort keys %ifile) {  # CHECK and LINK
    if (-r $f) {
      mylink("Inputs: ",$f,$ifile{$f}) ;
    } else {
      print "Missing Input $f !!!\n";
      die "ERROR: Missing $f (or possibly wrong Chem$Chem)\n" 
        unless( $f =~ /special/ or $f =~ /natso2/ );#natso2 not needed
    }
  }

  EMEP::Sr::generate_updated_femis(@$scenflag) if ($SR);
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

  my ($startdate,$enddate)=("$year-$mm1-$dd1","$year-$mm2-$dd2");
     $enddate=date2str($startdate." 1 day ago","%F") unless $dd2;
     ($startdate,$enddate)=("$CWFDATE[1]","$CWFDATE[2]") if $CWF;
# check if first/last met file exists
  foreach my $f ("$startdate","$enddate") {
    $f=date2str($f,$METformat);
    die "METFILE not found:\n\t$f\n" unless -e $f;
  }

# namelist with input parameters (to be read by Unimod.f90 and other modules)
  my $nml="";
  my %h=(%inml,%BENCHMARK);
  if(%BENCHMARK){
    # read nml template file
    $nml=EMEP::Sr::slurp("$ProgDir/config_BM-$GRID.nml");
    # fill in variables on the template file with corresponding $hash{key}
    %h=(%h,'runlabel1'=>"$runlabel1",'runlabel2'=>"$runlabel2",'METformat'=>"$METformat");
  } else {
    $nml="&INPUT_PARA\n"
        ."  iyr_trend = $iyr_trend,\n"
        ."  runlabel1 = '$runlabel1',\n"
        ."  runlabel2 = '$runlabel2',\n"
        ."  startdate = ".date2str($startdate ,"%Y,%m,%d,000000,\n")
        ."  enddate   = ".date2str($enddate   ,"%Y,%m,%d,000000,\n")
        ."  meteo     = '$METformat',\n"
        ."&end\n";
    # NML namelist options.
    foreach my $f ("config_$exp_name.nml","config_Outputs_$outputs.nml") {
      $nml.=EMEP::Sr::slurp("$ProgDir/$f");
    }
    # fill in variables on the template file with corresponding $hash{key}
    %h=(%h,'outdate'=>date2str($CWFDUMP[0],"%Y,%m,%d,000000,")
                     .date2str($CWFDUMP[1],"%Y,%m,%d,000000")) if $CWF;
  }
  $nml =~ s/(\s*\!.*|\s+$)//g;  # remove comments, tailing spaces
  $nml =~ s/\s*\n+/\n/g;        # & empty lines
  $nml.="\n";                   # restore newline at the end of last namelist
  # fill in variables on the template file with corresponding $hash{key}
  foreach my $k (keys %h) { $nml=~s:\$$k:$h{$k}:g; }
  # list mode setup variables
  $nml.="!". "-"x22 ." Model set-up ". "-"x22 ."\n";
  %h=('testv'=>"$testv",'Chem'=>"$Chem",'exp_name'=>"$exp_name",
      'outputs'=>"$outputs",'GRID'=>"$GRID",'MAKEMODE'=>"$MAKEMODE");
  foreach my $k (sort keys %h) { $nml.=sprintf "!  %-22s = '%s',\n",$k,$h{$k}; }
  # list of linked files in nml format for compatibility with future nml-only versions
  $nml.="!". "-"x22 ." Linked files ". "-"x22 ."\n";
  foreach my $k (sort keys %ifile) { $nml.=sprintf "!  %-22s = '%s',\n",$ifile{$k},$k; }
  open(TMP,">config_emep.nml");
  print TMP "$nml";
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

# Special CWF input files
  my $CWFINPUT="CWF? $CWF";
     $CWFINPUT.="\n".
     "  IC: ".(-e $cwfic?$cwfic:"Not found")."\n".
     "  BC: ".(-e $cwfbc?$cwfbc:"Not found")."\n".
     "  PL: ".(-e $cwfpl?$cwfpl:"Not found") if $CWF;

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
SR?  $SR
$CWFINPUT
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
     $last_sondes = "sondes_$year.csv" if ($CWF);
  print "LOOKING FOR LAST SITES $last_sondes\n";
  if ( -r $last_sondes ) {
    print "FOUND LAST sondes $last_sondes\n";
    system("tar cvzf $runlabel1.sites.tgz  sites.*");
    system("tar cvzf $runlabel1.sondes.tgz sondes.*");
  }

  if ($CWF) {
    my $old="EMEP_OUT.nc";
       $old=date2str($CWFBASE." 1 day","EMEP_OUT_%Y%m%d.nc") unless (-e "$old");; # 1st dump/nest
       $old=date2str($CWFBASE." 2 day","EMEP_OUT_%Y%m%d.nc") unless (-e "$old");; # 2nd dump/nest
    my $new=date2str($CWFBASE,$CWFIC);      # today's dump
    system("mkdir -p `dirname $new`; mv $old $new") if (-e "$old");
    if ($SR) {
      system("rm $cwfic") if (-e $cwfic);   # yesterday's dump
      $old="modelrun.finished";
      foreach my $task ('PBS_ARRAY_INDEX', 'PBS_ARRAYID', 'TASK_ID') {
        $new="../CWF_$CWFBASE/runsr_$ENV{$task}.finished" if $ENV{$task};
        system("mkdir -p `dirname $new`;echo $scenario >> $new")
          if (-e $old) && ($ENV{$task});
      }
    }
    # Pollen
    $old="POLLEN_OUT.nc";
    $old=date2str($CWFBASE." 1 day","POLLEN_OUT_%Y%m%d.nc") unless (-e $old);; # 1st dump/nest
    $old=date2str($CWFBASE." 2 day","POLLEN_OUT_%Y%m%d.nc") unless (-e $old);; # 2nd dump/nest
    $new=date2str($CWFBASE,$CWFPL);        # today's pollen dump
    system("mkdir -p `dirname $new`; mv $old $new") if (-e "$old");
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
  # links files from the original olcation (old) to
  # the new location (new) - generally the working directory.
  # Keeps track of all such linked files in list_of_files.
  my ($text,$old,$new) = ($_[0],$_[1],$_[2]);
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
  my ($cc,$pollut,$redn) = @_;
  # Initialise to 1.0:
  my($sox,$nox,$voc,$nh3,$testp,$co,$pm25,$pmco) = ("1.0")x8 ;
  given($pollut){
    when("AV"  ){ $voc = $nh3 = $redn };
    when("A"   ){ $nh3 = $redn };
    when("V"   ){ $voc = $redn };
    when("S"   ){ $sox = $redn };
    when("N"   ){ $nox = $redn };
    when("P"   ){ $pm25 = $pmco = $redn };
    when("NP"  ){ $nox = $pm25 = $pmco = $redn };
    when("SNP" ){ $sox = $nox = $pm25 = $pmco =  $redn };
    when("AN"  ){ $nh3 = $nox = $redn };
    when("SNAV"){ $sox = $nox = $nh3 = $voc = $redn };
   #when("BASE") then no change!
  }

  my $femisdat = slurp("$DataDir/femis.dat");
  die "ERROR!! No country Num for $cc!\n" unless defined(my $ccnum = $country_nums{$cc});

  # using 0 here as long as emissions are guaranteed to contain either
  # only anthropogenic or only natural emissions perl 'country'
  my $ss = 0; # 100 = antropogenic sectors (1..10)
              # 0 = all sectors
  $femisdat .= "$ccnum $ss  $sox $nox $voc $nh3 $testp $co $pm25 $pmco\n";
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
  }
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
