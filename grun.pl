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

require "flush.pl";

# <---------- "Pattern files - different on huge and gridur ---->
# These are changed by the  script accroding to the domain and
# number of processors selected.
# (but note that this script file is onlt for gridur these days)
#huge: @pattern_files = ( Makefile.pat Par_ml.pat );

@pattern_files = qw ( Par_ml.pat );

# <---------- start of user-changeable section ----------------->

#  --- Here, the main changeable parameters are given. The variables 
#      are explained below, and derived variables set later.-

$year = "1997";
( $yy = $year ) =~ s/\d\d//; #  TMP - just to keep emission right
# NEW: iyr_trend can be set to meteorology year or arbitrary year, say 2050
$iyr_trend = $year;  
#$iyr_trend = "2010";     # For CLE runs, SR runs, etc.

print "Year is $yy YEAR $year\n";

if ( $year == 2000 ) {
  $MetDir = "/work/mifapw/metdata/$year" ;
} elsif ( $year == 1997 ) {
  $MetDir = "/work/mifapw/metdata/$year" ;
} elsif ( $year == 1995 ) {
  $MetDir = "/work/mifaab/metdata/$year" ;
} elsif ( $year == 1990 ) {
  $MetDir = "/work/mifads/metdata/$year" ;  # Needed for 1997
} else {
  $MetDir = "/work/mifahf/metdata/$year" ;
}

#---  User-specific directories (changeable)

$DAVE        = "/home/u2/mifads";      
$JOFFEN      = "/home/u2/mifajej";      
$HILDE       = "/home/u5/mifahf";      
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
# Also, list emissions to be used:

if ( $OZONE ) {
    #$OZONEDIR    = "$HILDE/BC_data/CTM2_data/50Data"; 
     $OZONEDIR    = "$HILDE/BC_data/LOGAN_O3_DATA/50Data_900mbar"; 
    #$OZONEDIR    = "$HILDE/BC_data/Fortuin_data/50Data"; 
     @emislist = qw ( sox nox nh3 co voc pm25 pmco ); 
     $testv       = "rv1_9_5";
     $runlabel1    = "CLE";   # NO SPACES! SHORT name (used in CDF names)
     $runlabel2    = "${testv}_CLE_$year";   # NO SPACES! LONG (written into CDF files)

} elsif ( $ACID ) {
     $OZONEDIR    = "$HILDE/BC_data/EMEPO3_rv147";
     @emislist = qw ( sox nox nh3 pm25 pmco ) ;
     $testv       = "rv1_7_3";
} 
#$H2O2DIR     = "$HILDE/BC_data/EMEPH2O2_rv147";     # Needed for both acid and ozone
$H2O2DIR     = "$HILDE/BC_data/EMEPH2O2_rv1.5.1oxlim";# Needed for both acid and ozone
$version     = "Unimod" ;  
$subv        = "$testv" ;                  # sub-version (to track changes)
$Case        = "DSTEST" ;                   #  -- Scenario label for MACH - DS
$ProgDir     = "$USER/Unify/$testv"; # input of source-code
$MyDataDir   = "$USER/Unify/MyData";          # for each user's femis, etc.
$DataDir     = "$DAVE/Unify/Data";      # common files, e.g. ukdep_biomass.dat
$PROGRAM     = "$ProgDir/$version";         # programme
$WORKDIR     = "$WORK{$USER}/Unimod.$testv.$year";    # working directory
$femis       = "$MyDataDir/femis.dat";      # emission control file

# Fixed
#$emisdir     = "$SVETLANA/Unify/MyData/emission";   # emissions directory
#Latest: (not used for PM though).
#$emisdir     = "$HILDE/emis/trends2003"; # emissions directory
#$emisyear    = "$emisdir/emis${yy}-V3";    # emissions
$emisdir     = "$PETER/Unify/MyData/"; # emissions directory
$emisyear    = "$emisdir/2010_CLE_2000";    # emissions
$timeseries  = "$DAVE/Unify/D_timeseries";   # New timeseries (ds 14/1/2003) 

# Specify small domain if required. 
#                 x0   x1  y0   y1
@largedomain = (   1, 170,  1, 133 ) ;
#@smalldomain = ( 101, 140, 51,  90 ) ;      # (changeable)
#@smalldomain = (  71, 150, 31, 100 ) ;      # (changeable)
#@smalldomain = (  95, 115, 46, 66 ) ;      # ERROR search (changeable)
#@smalldomain = (  36, 160, 11, 123 ) ;      # (changeable)
#@smalldomain = (  20, 167,  1, 122 ) ;    # OSPAR/HELCOM domain
@smalldomain = (  36, 167, 12, 122 ) ;    # EMEP domain
#@smalldomain = (  36, 130, 31, 123 ) ;      # (changeable)
#@smalldomain = (  39, 120, 31, 123 ) ;      # (changeable)
#@smalldomain = @largedomain ;     # If you want to run for the whole domain, 
                                    # simply uncomment this 

$RESET        = 0   ;  # usually 0 (false) is ok, but set to 1 for full restart
$COMPILE_ONLY = 0   ;  # usually 0 (false) is ok, but set to 1 for compile-only
$INTERACTIVE  = 0   ;  # usually 0 (false), but set to 1 to make program stop
                       # just before execution - so code can be run interactivel.

$NDX   =  8;           # Processors in x-direction
$NDY   =  4;           # Processors in y-direction
if ( $INTERACTIVE ) { $NDX = $NDY = 1 };


@month_days   = (0,31,28,31,30,31,30,31,31,30,31,30,31);
$month_days[2] += leap_year($year);

$mm1   =  1;       # first month
$mm2   = 12;       # last month
$NTERM_CALC =  calc_nterm($mm1,$mm2);

$NTERM =   $NTERM_CALC;    # sets NTERM for whole time-period
  # -- or --
# $NTERM =  12;       # for testing, simply reset here

  print "NTERM_CALC = $NTERM_CALC, Used NTERM = $NTERM\n";


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


#--- Change to WORKDIR
chdir "$WORKDIR"; 

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
    # Ozone and H2O2 data are needed for all model versions:
    # (Assumes all ozone data have same name, i.e. ozone.mm)

    $old = sprintf "$OZONEDIR/ozone.%02d", $mm ;
    $new = sprintf "ozone.%02d", $mm ;
    mylink( "Linking BCs:", $old,$new ) ;

    $old = sprintf "$H2O2DIR/h2o2.%02d", $mm ;
    $new = sprintf "h2o2.%02d", $mm ;
    mylink( "Linking H2O2 BCs:", $old,$new ) ;

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

    $old   = "$MyDataDir/femis.dat" ;
    $new   = "femis.dat";
    mylink( "Femis  ", $old,$new ) ;

    $old   = "$MyDataDir/splits/pm25split.defaults.$Case" ;
    $new   = "pm25split.defaults";
    mylink( "Split pm25", $old,$new ) ;

    $old   = "$MyDataDir/splits/vocsplit.defaults.$Case" ;
    $new   = "vocsplit.defaults";
    mylink( "Split voc", $old,$new ) ;

    $old   = "$MyDataDir/splits/vocsplit.special.$Case" ;
    $new   = "vocsplit.special";
    mylink( "Split voc", $old,$new ) ;

foreach $poll  ( @emislist  ) {
#  if ( $poll =~ /pm/ ) {  # CRUDE FIX FOR NOW
#   $old   = "$SVETLANA/Unify/MyData/emission/em2000/grid$gridmap{$poll}" ;
#  } else {
   $old   = "$emisyear/grid$gridmap{$poll}" ;
#  }
   $new   = "emislist.$poll";
   mylink( "Emis $poll : ", $old,$new ) ;

   #rv141:$old   = "$emisdir/Monthlyfac.$poll" ;
   $old   = "$timeseries/MonthlyFac.$poll" ;
   $new   = "MonthlyFac.$poll";
   mylink( "MonthlFac ", $old,$new ) ;

   $old   = "$timeseries/DailyFac.$poll" ;
   $new   = "DailyFac.$poll";
   mylink( "DailyFac ", $old,$new ) ;
} 

# Surface measurement sites
# From $DAVE directory for these
    $old   = "$DataDir/sites.rep03" ;
    $new   =  "sites.dat";
    mylink("Sites ",  $old,$new ) ;

# Sondes
    $old   = "$DataDir/sondes.rep03" ;
    $new   = "sondes.dat";
    mylink( "Sondes", $old,$new ) ;

# Forest data
    #ds $old   = "$o3dir/forest.pcnt" ;
    #ds $new   = "forest.pcnt";
    $old   = "$DataDir/forests.tf2" ;
    $new   = "forest.tf2";
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

    $old   = "$DataDir/landuse.nov2003" ;  #ds rv1_9_4 change
    $new   = "landuse.dat";                #ds rv1_9_4 change
    mylink( "Landuse ", $old,$new ) ;

 # TMP LOCATION for some datafiles : MyDataDir
#dsforeach $datafile ( qw ( Volcanoes.dat ukdep_gfac1.dat ukdep_gfac2.dat ukdep_biomass.dat ) ) {
foreach $datafile ( qw ( Volcanoes.dat tf2_gfac1.dat tf2_gfac2.dat tf2_biomass.dat ) ) {
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

foreach $exclu ( @exclus) {
    print "starting $PROGRAM with 
        NTERM $NTERM\nNASS $NASS\nEXCLU $exclu\nNDX $NDX\nNDY $NDY\nIYR_TREND $iyr_trend\nLABEL1 $runlabel1\nLABEL2 $runlabel2";

if ( $INTERACTIVE ) {
  die " -- INTERACTIVE: can now run for inputs: NTERM, NASS,  EXCLU, NDX, NDY 
  with mpirun -np 1 prog.exe << XXX
  $NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\nXXX"; 
}

    #open (PROG, "|mpirun -np $NPROC $PROGRAM") || 
    open (PROG, "|mpirun -np $NPROC $LPROG") || 
               die "Unable to execute $LPROG. Exiting.\\n" ;

    # open (PROG, "|mp
    #   open (PROG, "|$PROGRAM") || 
    #          die "Unable to execute $PROGRAM. Exiting.\\n";
    
    #print PROG "$NTERM\n\n$NASS\n$exclu\n$NDX\n$NDY\n\iyr_trend\n";
    #print PROG "$NTERM\n$NASS\n$exclu\n$NDX\n$NDY\n\iyr_trend\n";
    print PROG "$NTERM\n$NASS\n$iyr_trend\n${runlabel1}\n${runlabel2}\n";
    close(PROG);
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

