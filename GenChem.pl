#!/usr/bin/env perl
use strict;
use warnings;
# =====================================================================
#  GenChem - a perl script to read in equations and output
#  production and loss terms for fortran programmes.
#
# Dave Simpson    1998-2013:
# Includes: species checks, atom-balance,  use of shorthands,
#  rate-simplication, various user-changeable options.
#
# Modified by Garry Hayman and Dave, 2009 for new-style input files
# and EMEP rv3_3 model version.
#  =====================================================================
#
#  This  script starts with a section of "user-definable" text-strings
#  which can be used to write different output formats, e.g. for Fortran
#  90, C, 1-D or 5-D models. The main programme follows and is very short,
#  simply calling the 6 subroutines that do the work.
#
# >>>>>>>>>>>>>>>>>   Definitions : model specific  <<<<<<<<<<<<<<<<<<<
#                     May be changed by user
#
# Define the text used at the start and end of each production or loss
# term. This is defined here for the EMEP 3D model, and probably
# should be changed for other models
#
our $DEBUG = 0;
my $DIMFLAG= "k";       # could be ix,iy,iz for 3-D models if wanted
my $TEMP   = "XT";      # special variable
my $Krange = "KCHEMTOP:KMAX_MID";

my $PROD = "P" ;           # Could be P, P(iq), P(i), whatever
my $LOSS = "L" ;
my $CPNEW = "xnew";
my $CPDIM = "";   # could be ",k" in some routines

my $LOOP_HEADTEXT =  "";    #  Could be "DO\n" -loop
my $LOOP_TAILTEXT =  "";    #  e.g.     "END DO\n";

my $CHEMEQN_TEXT =  "
      xnew(SPEC)=  ( xold(SPEC) + dt2 * P) /(1.0 + dt2*L )\n";
      #if( xnew(SPEC) < 0.0 ) print *, \"NEGNEG!!! SPEC\", xold(SPEC), P, L\n";
      #xnew(SPEC)=  ( xold(SPEC) + dt2 * P) /(1.0 + dt2*L )\n";

my $CHEMEQN_NOLOSS  =  "
      xnew(SPEC)=  xold(SPEC) + dt2 * P\n";

my $CHEMEQN_NOPROD  =  "
      xnew(SPEC)=  xold(SPEC) / ( 1.0 + dt2 * L )\n";

my $CHEMEQN_NOTERMS =  "!Nothing to do for SPEC! xnew(SPEC)=  max(0.0, xold(SPEC))\n";

# From the original towstep...
#
# $CHEMEQN_TEXT =  "      X(iq,SPEC,2)=
#      >     (4.0*X(iq,SPEC,1)-X(iq,SPEC,0)+2.0*Ldt*P)/
#      >     (3.0 + 2.0*Ldt*L )\n";

open(LOG,">Log.GenOut");    # An extensive log - should be read!

# flag to cope with new-style files from Garry. New files allow RO2_POOL to be used
# And grouped style from Dave, added Mar 2010
my $old_style = 0;
our $grp_style = 0;
our %grp;   # Will be hash of arrays,  e.g. $grp{"NOX"} = [ "NO", "NO2" ]

our @emis_specs = ();
our @emis_files = ();
our @bionat_specs = (); # For isoprene, terpenes, soil-NO etc.

# Some descriptive variables for later output
my ( $line, $linenum, $rate ) ;  # Input line
my ( $shl, $adv, $svol, $tot, $tracer )   = ( 0, 1, 2, 3, 4 ) ;
my $SlowSpec = 9999;   # Index of slowly-reacting species
my ( %nspec_txt   )  = ( $adv => "NSPEC_ADV",   $shl => "NSPEC_SHL" );
our ( @species_tot ); # Just used in checking drydep now

#2010: reads groups and headers from GenIn_species.csv
my( %groups, @headers );
my @ro2pool = (); # used with new-style GenIn.species from Garry

my ( $naerosol, $first_semivol, $last_semivol )   = ( 0, 0, 0 ) ;
my ( @species, @nspecies, @speciesmap ) = ();     # e.g. species[TOT][n]
my ( %extinc, @CiStar, @DeltaH ); # Aerosol params

# Arrays for left and righthand side of equations:
my (@lhs,  @rhs) = ();               # Initialise lhs and rhs species
my (@lhtracer,  @rhtracer) = ();     # Initialise lhs and rhs tracers
my (%lhfactor,  %rhfactor ) = ();    # Initialise lhs and rhs factors
my (%prod,  %loss ) = ();            # Initialise
my ($neqn, $nrct, $nrctroe, $nrcfunc ) = (0,0,0,0);  # rct = Temp-dependant
my ( $nddep,  $nwdep ) = (0,0);
my ( $ddep_txt, $wdep_txt );
#QUERY? my ( $ddep_txt, $wdep_txt ) = ("","");
my ( @rate_label, @rct, @rcttext, @rctroe, @rctroetext,
  @rcmisc, @rcmisctext ) = ();
my ( %atoms, %molwt, %count, %nmhc, %shorthand ) = ("", "", "", "" );

# Common statements for output:
my $CONTLINE = "  &\n     ";      # Fortran'90
my $INTPUB   = " integer, public, parameter :: " ;
my $LOGPUB   = " logical, public, parameter :: " ;
my $HLINE = "!-----------------------------------------------------------\n";
my $TLINE = "!>_________________________________________________________<\n";

#======================================================================
# CONSIDER REMOVING!
my ( $rct )    = qw ( rct );
my %description = ($adv => "Advected species",
                   $shl => "Short-lived (non-advected) species ",
                   $tot => "All reacting species ",
                   $rct => "Rate-coefficients - temperature dependant");

#  These rates require different used variables. The ChemFunctions_ml
#  is only called if the more complex functions are needed.
#  but for simplicity we define it anyway.... (ESX check)
my $txtrcfunc = "";
$txtrcfunc = "\n  use ChemFunctions_ml       ! => kaero, RiemerN2O5\n" ; #ESX if $nrcfunc>0;
#printall ("TXTRCFUNC needed? $nrcfunc   TXT:$txtrcfunc");
my %UsedVariables = (
   $rct    => " $txtrcfunc
  use AeroFunctions     ! => UpdakeRate, cMolSpeed
  use Setup_1dfields_ml ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use ChemSpecs_tot_ml  ! => PINALD, .... for FgasJ08
  use ModelConstants_ml, only: KMAX_MID,KCHEMTOP,DebugCell,DEBUG,AERO");

#======================================================================

read_species();    # Assigns to $species[nspec] array, also calculates
                   # $molwt{$spec} and $count[$spec][$atom].
print "DONE READ_SPECIES\n";

print_species();
print "DONE   PR_SPECIES\n";

read_shorthand();  # Reads and expands shorthand, e.g. KHO2RO2
print "DONE READ_SHORTS\n";

read_reactions();  # Reads reactions, does carbon balance etc.
print "DONE READ_REACTIONS\n";

output_prod_loss();  # Output production and loss file : GenOut_reactions.inc
print "DONE OUT_PROD\n";

print_rates($rct, $nrct,\@rct,\@rcttext) ;

printmap_dep();
printmap_ext();
print "QQQQQQQQQQQQQQQQQQQQQQQQQQ ; \n";
foreach my $gg ( keys %grp ) {
  print "  TESTQQ $gg group is: @{ $grp{$gg} }\n";
}
print_groups();   # DDEP_OXNGROUP etc.
print_emisstuff("File",@emis_files );     #  "nox ", etc.
print_emisstuff("Specs",@emis_specs );     #  "nox ", etc.
print_emisstuff("BioNat",@bionat_specs );     #  "nox ", etc.
print_femis(@emis_files );     #  "nox ", etc.
print_emislist(@emis_files );     #  "nox ", etc.

#print_speciesmap() ;


# =====================================================================
#######################################################################
sub read_species {
  print LOG "............... read_species ..........................
          reads from GenIn.species, assigns nspec, and species[nspec],
          and calls count-atoms to calculate weights and atom counts.
          Note we have to be careful with case here. The model name
          for the species is in uppercase, but the chemical forumla
          read in needs correct case, e.g. Na, Cl, etc.
          New:
           The input typ can be:       0 for non-advected species
                                       1 for advected species
                                       2 for gas-particle (soa) species
  ...................   \n";

  open(F,"GenIn.species") or die "File GenIn.species not found\n";

  while ( $line = <F> ) {       # Reads a line
    chomp($line) ;         # Remove newline character
    next if $line =~ /^\*/ ;   # Skip comments
    next if $line =~ /^\"\*/ ;   # Skip comments
    next if $line =~ /^\s*$/ ; # skip empty lines
    if ( $line =~ /^Spec/ ) {
      @headers = split(/,/,$line);
      print join(" ",@headers);
      next;
    }

    $SlowSpec = $nspecies[$tot] + 1 if $line =~ /^\"#SLOW/ ;
    next if $line =~ /^"#/ ; # would be OLD_STYLE of SLOW

    my ( $spec, $typ, $formula, $comment, $in_ncarbon, $class );
    my ( $extinc, $cstar, $DeltaH ) = 0.0x3; # Aerosol params
    my ( $groups, @groups );
    @groups = ();
    my $in_rmm = "-";   # init to "-" for old_style
    $line =~ s/xx/-/g; # Oocalc didn't like -
    my @r = split(/,/,$line);
    $spec     = $r[0];
    $typ      = $r[1];
    $formula  = $r[2];
    $in_rmm   = $r[3];
    my $dry   = $r[4];
    my $wet   = $r[5];
    $extinc   = $r[6];
    $cstar    = $r[7];
    $DeltaH   = $r[8];
    $groups   = $r[10];
    my $tgroups  = $r[ find_rec("Groups",@headers ) ];
    #die "TESTED T$tgroups should be ten\n";

    # Assign species to groups if given:
    unless ( $groups =~ "-" ) {
      print "PG SPEC $spec G $groups\n";
      process_groups($spec, $groups, $dry, $wet );
      #print "PPPPPPPPPPPPPPPPPPPPPPPPPP ; \n";
      #foreach my $gg ( keys %grp ) {
      #print "  TESTPP $gg group is: @{ $grp{$gg} }\n";
      #}
    }
    process_alldep ( $dry, $wet, $spec ) unless $dry eq "-" ;

    #???? push(@ro2pool,$spec) if $class eq "peroxy";
    $spec  =~ tr/a-z/A-Z/;   # Upper case

    $nspecies[$tot]++;
    my $n = $nspecies[$tot];
    $species[$tot][$n ]    = $spec ;
    $species_tot[$n]       = $spec ;
    $speciesmap[$tot][$n]  = $n ;
    $extinc{$spec} = $extinc unless ($extinc =~ /0|-|^$/);
    $CiStar[$n]  = $cstar;
    $DeltaH[$n]  = $DeltaH;
    printall("SPECIES READ => $n $spec N $n DH $DeltaH; ");

    if ( $typ == 2 ) {
      $naerosol++ ;
      if ( $naerosol == 1 ) { $first_semivol = $n };
      $last_semivol = $n ;
      $typ = 1;  # Reset for advection below
    }

    # also for advected, short-lived, semi-vol types:

    $nspecies[$typ]++;
    $n = $nspecies[$typ];
    $species[$typ][$n ]    = $spec ;
    $speciesmap[$typ][$n]  = $n ;

    print "LINE $line TYPE $typ N $n\n";
    if (is_integer($formula) or is_float($formula)) {
      $molwt{$spec} = $formula;
      print LOG "INPUT MOLWT: $spec= $molwt{$spec}\n";
    } else { 
    # Get $molwt and $count[$spec][$atom]
      count_atoms($spec,$formula)
    }

    if ($molwt{$spec} == 1 and $in_rmm ne "-" ) {
      $molwt{$spec} = $in_rmm
    }

    if (is_float($in_rmm) or is_integer($in_rmm)) {
      $molwt{$spec} = $in_rmm;
      print LOG "INPUT MOLWT: $spec= $molwt{$spec}\n";
    }

    unless (defined($count{$spec}{"C"})) {
      if ($count{$spec}{"C"} == 0 and $in_ncarbon != 0)
        { $count{$spec}{"C"} = $in_ncarbon }
    }
    #AMVB: NMVOC group
    push @{$grp{'NMVOC'}},uc($spec) 
      if($typ eq $adv)and($count{$spec}{"C"}>0)and
        (uc($spec)ne"CO")and(uc($spec)ne"CH4");

  }
  close(F);


  print "PPPPPPPPPPPPPPPPPPPPPPPPPP ; \n";
  foreach my $gg ( keys %grp ) {
    print "  TESTPP $gg group is: @{ $grp{$gg} }\n";
  }
} # end of sub read_species

########################################################################
sub read_shorthand {
print LOG "\n\n............... read_shorthand ........................
         reads from GenIn.shorthand, e.g.:
         XT           temp(iq)
         H2O          H2O(iq)
         FH2O         (1.+1.4e-21*h2o*exp(2200./XT))
         and creates shorthand with expansions, e.g. \$shorthand{XT}= temp(iq)
.........................................................................  ";
        my ( $line,$pattern,$expanded);
#       global: establishes $shorthand{$pattern}

        open(F,"<GenIn.shorthand");
	printall ("\n\n\nThe following SHORTHANDS were processed: \n\n");

        while ( $line = uc <F> ) {         # Reads a line, make upper-case
                chomp($line) ;              # Remove newline character
                next if $line =~ /^\*/ ;   # comment line (starts with *): skip
                next if $line =~ /^\s*$/ ; # empty line : skip

                my @input = split(/\s+/,$line,3);   # split line into 3 columns
                                                 # (3th column contains any comments)
		$pattern  = $input[0];
		$expanded = $input[1];

		# Expand any shorthand in $expanded - eg. XT to temp(iq)
		$expanded = expand_shorthands( $expanded );
		#foreach my $shorthand (keys(%shorthand)){
		#   if ($expanded  =~ /\b $shorthand \b/x ) {
		#	$expanded  =~ s/$shorthand/$shorthand{$shorthand}/g;
		#   }
		#}

		$shorthand{$pattern} = $expanded ;  # --> e.g. $shorthand{HO2RO2}=1.0e-11
		printall ("$pattern  - >   $shorthand{$pattern}\n") ;
        }
        close(F), printall ("\n\n ...................     \n") ;
} # end of sub read_shorthand

#########################################################################
sub read_reactions {
print LOG "\n\n............... read_reactions .........................
         reads from GenIn.reactions, assign nspec, and species[nspec]
........................................................................\n" ;
         my ($i, $term);

	 $linenum = 0;     # Line No. of input file, for checking
         open(F,"GenIn.reactions") or die "FAILURE: open of GenIn.reactions\n";

         while ( $line = uc <F> ) { # Reads a line at a time from F, upper-case
                chomp($line) ;       # Remove newline character
		$linenum++ ;

        	next if $line =~ /^\*/ ;        # comment line (starts with *)
        	next if $line =~ /^\s*$/ ;      # blank line (skip)
        	next if $line =~ /YA\(\d+\)/ ;  #  Garry's SOA
        	$line =~ s/YG\(\d+\)// ;        #  Garry's gas-phase, skip yg() terms
		print " PRE $line \n";
        	$line =~ s/\;.*/;/ ;             #  Remove comments
		print " POST $line \n";
	        if( $line =~ /emisfile/i ) {
			printall( "Process EMISFILES??: emisfile $line \n" );
 			process_emisfile( $line ) ; # rcemis:NO2
			#printall( "Process??: K $k \n" );
			next;  # Not a reaction
		 }

		$neqn++ ;
	        printall ("."x66 . "\n");
	        printall (" "x19 . "Processing Line $linenum Equation $neqn :\n $line\n\n");
		# check for improper equations (e.g. forgotten comment lines without *)
		die "Not a proper equation! No = sign! " unless $line =~ /=/ ;

		my @terms = split ' ',$line ;      # split $line into terms

		$rate = $terms[0];


		# Split rates into constants, variables, emis terms, etc:

		$rate_label[$neqn] = define_rates() ;

		# Process LHS (reactants)

	 	$i = 0 ;
		@lhs = (); @lhtracer = (); %lhfactor = ();
		#while( ( $term = $terms[++$i]) ne "=" ) {   # Loop through $terms
			#	printall("TESTING TERM $i $term \n" );
			#}
		#$i = 0 ;
		while( ($term = $terms[++$i]) ne "=" ) {   # Loop through $terms
			printall("TESTING TERM $i $term \n" );
			if ( $term  eq "+" ){
				;
			} elsif ( $term =~ /\+/ ){     # Check for easy-to-make error
				die "ERROR!! Found + sign without surrounding space. Correct!\n";

			} elsif ( $term =~ /{.+}/ ){     # Tracers specified with {XX}
				push(@lhtracer,$term) ;
				count_tracer($term);    # Count atoms

			} elsif ( $term =~ /\[(.+)\]/ ){     # Tracers specified with {XX}
				my $mult = $1; # matches contents
				push(@lhtracer,$term) ;
				count_tracer($term);    # Count atoms
				$rate_label[$neqn] .= "*$CPNEW($mult)";
				print "TRACER MULT $term => $rate_label[$neqn]\n";
			} else {
				print "LHS PUSH  $term => @lhs \n";
			    	push(@lhs,$term);      # Add to lhs array
			}
		}
	 # process_lhs sets losses fills in loss matrix, and also returns
	 # $flux = production flux from left hand side...

		my $flux = process_lhs() ;

		# Process RHS (products), this time one term at a time

		my $factor = 1;
		@rhs = (); @rhtracer = (); %rhfactor = ();
		while ( ($term=$terms[++$i]) ne ";" ){
			if ( &is_float($term)  ){  # matches any number (inc. eE type, but not dD)
			  $factor = $term ;
			} elsif ( $term eq "+" ){
			   ;
			} elsif ($term =~ /{.+}/ ){ # Matches {} - tracers
			   push(@rhtracer,$term) ;  # Add tracer
			   count_tracer($term);    # Count atom
			   if( $factor != 1) {printall("Used RH factor $factor for $term\n")};
			   $rhfactor{$term} = $factor, $factor = 1;
			} else{
			   die "Error:
			    do not include multiplication signs on RHS
			    :: $line\n" if( $term eq "*" );
			   push(@rhs,$term);      # Add to rhs array
	printall("CHECK_MATCH FROM RHS  RR${term}RR ");
			   check_match( $term ) ;  # Check that species exists
			   process_rhs( $term,$factor,$flux ) ;
			   if( $factor != 1) {printall("Used RH factor $factor for $term\n")};
			   $rhfactor{$term} = $factor, $factor = 1;
			}
		}
		carbon_count_etc() ;
	}
 	close(F);
} # end of sub read_reaction

########################################################################
sub count_atoms {
##      Counts the number of each atom of carbon, nitrogen, etc.,
##      - should one day add the whole periodic table to catch all
##        possibilities ! Note that for species whose weight could
##        not be determined a value of 1 is used.
#######################################################################
my($spec, $formula) = @_;
my($num) ;
# Global - defines %atoms, calculates $count[$spec][$atom]


%atoms = ("C" => 12, "H" => 1, "N" => 14, "O" => 16, "S" => 32 );

# Perl note - the following uses a match with parentheses =~ /(x)(y)/
# perl sets $1 = x and $2 = y in this case, and seems to reset $2 to ""
# if y doesn't exist. Hmm... useful, but I don't know the rules used.

	#print LOG "Count_atom: $spec = $formula \t\t ";
        print LOG "Atoms: \t ";
        # look for H12, C3, etc., first, count  the atoms.
        foreach my $atom ( keys(%atoms) ){
                $num = 0;
                while ( $formula =~ /($atom)([\d]*)   # Looks for e.g. N, N2
                                    (?![a-z])         # Not followed by a-z, (eg. Na)
                                    /gx ){            # g: through line, x:allow comments
                                                      # (like these!)
			# Was $2 == 0 ....
                        if($2 eq "" ) {    # If we have no number (e.g. N)
                                $num+= 1 ;
                        } else{            # we had some number! (e.g. N2)
                                $num+= $2 ;
                        }
                }
                $molwt{$spec} += ( $num *  $atoms{$atom} );
                $count{$spec}{$atom} = $num ;
                if ( $num ) {print LOG " $atom=$count{$spec}{$atom}"} ;
        } # atom loop

        # Test for NMHC, contains only C(+H), not O, S, N. Exclude CH4
        $nmhc{$spec} = 0 ;
        if ( $count{$spec}{"C"} >= 2.0 ) { $nmhc{$spec} = 1 };
        if ( $count{$spec}{"O"} >= 1.0 ) { $nmhc{$spec} = 0 };
        if ( $count{$spec}{"S"} >= 1.0 ) { $nmhc{$spec} = 0 };
        if ( $count{$spec}{"N"} >= 1.0 ) { $nmhc{$spec} = 0 };

	if ( $molwt{$spec} == 0 ) {
		print LOG "; NOT RECOGNISED ! " , $molwt{$spec}=1;
	}

        print LOG "\t Weight = $molwt{$spec}\n";  # End of scan
} # end of sub count_atoms ############################################

#########################################################################
sub define_rates {
#  ............... define_rates ..........................
#  splits rate-coefficients into proper constants, variables with
#  dimension DIMFLAG (=iq here, could be ix,iy if wanted), emission terms, etc.
#  Uses some cryptic perl notation, but it is effective!
#
#      Perl memos:
#       =~ is pattern matching operator, which uses following
#      \d = any digit   ([0-9])
#      + = one or more times, * = zero or more times, ? = zero or more
#          of one character only
#######################################################################
	my $k ;
	my $oldrate = $rate ;  # oldrate for Log

# First, we check if rate constant is one of the shorthands, e.g. KHO2RO2. If so
# we substitute the value (e.g. 1.0e-11)  into rate and carry on. To make
# the search easy, we first add spaces around operators

	printall("-------------------\n START Rate $oldrate\n");
	$rate =~ s/(\*|\+|\-\/)/ $1 /g ;

	$rate = expand_shorthands($rate);
	printall("Expand shorthands:  $rate \n");

# Then, we see what type of rate we have. We exclude exponentials
# and temperature(=XT) and see if any other text is present. If not
# we have either a pure number or a complex expression.

# For troe expressions, we replace the Fc term by its logarithm
# and replace the very common 300/xt by t300. Both of these
# are optimisations suggested by Steffen.

	my $troe  = 0 ;
	my $rctest = $rate;
	$rctest =~ s/EXP//g;
	$rctest =~ s/XT//g;
	$rctest =~ s/TINV//g;
	#ESX if ( $rctest =~ /[a-dfzA-DF-Z]/ ) { $misc = 1 } ;
	$troe = 1 if $rctest =~ /TROE/ ; #ESX ) { $troe = 1, $misc = 0 } ;
	if ( $troe ){ printall( "TRUETRO $troe for $rate \n")};

	$rate =~ s/\s+//g ; # Get rid of white space again
	# $rate =~ s/([\d.]+)[dD]([-+\d]+)/$1E$2/ ;
	printall("Type? $rate TROE $troe \n");

# Look for exponentials. Lower case them to keep F happy, and maybe later to
# help optimise.

	$rate =~ s/\b EXP \b/exp/xg ;              #  Change EXP to exp
	$rate =~ s/([\d\.]) ([E]) ([\+\-]?\d)/$1e$3/x ;     #  Change E to e
	$rate =~ s/([\d\.]) ([dD]) ([\+\-]?\d)/$1e$3/x ;     #  Change D/d to e
	$rate =~ s/ (\.)(\D)/\.0$2/gx ;    # Change e.g. 1.e14 to 1.0e14
	printall( "Clean: RATE $rate \n" );

# Now we continue to process, looking for pure numbers first:
    {   local $^W=0;   # Switch off -w flag for number checking

	if ( $rate =~ /[\*\+]/ )  {
                simplify_rate($rate) ;     # Simplify any arithmetic expressions:
	}

	if ( &is_float($rate) ) {
		$k = $rate ;
	} elsif ( $rate =~ /rcemis/i ) {
		printall( "Process??: RATE $rate \n" );
		$k = process_emis( $rate ) ; # rcemis:NO2
		printall( "Process??: K $k \n" );
	} elsif ( $rate =~ /rcbio/i ) {
		printall( "Process??: BIO RATE $rate \n" );
		$k = process_emis( $rate ) ; # rcbio:NO
		printall( "Process??: BIO K $k \n" );
	} elsif ( $rate =~ /aqrck/i ) {
		$k = $rate ;
	} elsif ( $rate =~ /DJ/i ) {   # e.g. DJ(iq,2)
		$rate =~ s/DJ/rcphot/;
		$k = $rate ;
#ESX	} elsif ( $rate =~ /DRY|WET/i ) {   # e.g. DJ(iq,2)
#ESX		$k = "$rate($DIMFLAG)" ;
#ESX
#ESX	} elsif ( $rate =~ /TROE/ or $rate =~ /KMT/ ) {   #  Troe expression - T & P
#ESX		$nrcmisc++ ;
#ESX		$k = "rcmisc($nrcmisc,$DIMFLAG)";
#ESX		$nrctroe++ ;
#ESX		$rate = lc($rate) ;
#ESX		printall( "A08TROE->RCMISC = $nrcmisc RATE $rate\n" );
#ESX		$rcmisc[$nrcmisc] = lc($rate) ;
#ESX		$rcmisctext[$nrcmisc] = "rcmisc($nrcmisc,:)" ;    # older: $k ;
#ESX	} elsif ( $rate =~ /(^[A-Z])/ ) {   # e.g. KRO2HO2,DJ_2  --> DJ_2(iq)
#M00		if ( $rate =~ /[\*\+\-]/ )  {
#M00                     die "Error in $rate, line $linenum :
#M00                     Rate-coefficients cannot begin with variable names, (e.g. K*3.0),
#M00                     put variable last (e.g. 3.0*K) or use GenIn.shorthands !\n" ;
#M00		}
#M00		$k = "$rate($DIMFLAG)"    ;
#ESX		$k = $rate  ;

	} elsif  ( $rate =~ /^_FUNC_/ ) {   # Complex expression
		#ESX $nrcmisc++ ;

		#ESX$k = "rcmisc($nrcmisc,$DIMFLAG)";

		$nrcfunc++; 
		$nrct++;  #ESX
		$rate =~ s/_FUNC_//;
		printall( "Function->RCMISC = $nrct RATE $rate\n" );
		#ESX$rcmisc[$nrcmisc] = lc($rate) ;
		#ESX$rcmisctext[$nrcmisc] = "rcmisc($nrcmisc,:)" ;    # older: $k ;

		$rct[$nrct] = $rate ;
		$rcttext[$nrct] = "rct($nrct,:)";    # older ,$DIMFLAG)";
		$k = "rct($nrct,$DIMFLAG)";

#ESX	} elsif  ( $misc == 1 ) {   # Complex expression
#ESX		$nrcmisc++ ;
#ESX		printall( "MRCMISC = $nrcmisc\n" );
#ESX		$k = "rcmisc($nrcmisc,$DIMFLAG)";
#ESX		$rcmisc[$nrcmisc] = lc($rate) ;
#ESX		$rcmisctext[$nrcmisc] = "rcmisc($nrcmisc,:)" ;    # older: $k ;
	} else{ # Variable coefficients --> e.g. rc(112,iq)
		# For optimum speed, check for previous_rates ........   ;
		for (my $i=1; $i<=$nrct;$i++){
		    if ( $rate eq $rct[$i] ){
			$k = "rct($i,$DIMFLAG)";
			return ($k);
		    }
		}
		# New coefficient needed...
		$nrct++;
                $nrctroe++ if $rate =~ /TROE/ ;  #  Troe expression - T & P
		$nrcfunc++ if ( $rate =~ /TROE/ or $rate =~ /KMT/ ) ;
                printall("ADDING TROE $rate") if $rate =~ /TROE/ ;  #  Troe expression - T & P
                printall("ADDING FUNC $rate") if $rate =~ /TROE/ or $rate =~ /KMT/ ;  #  Troe expression - T & P
		# if ( $rate =~ /[\*\+\-]/ )  { simplify_rate($rate) };
		$rct[$nrct] = $rate ;
		$rcttext[$nrct] = "rct($nrct,:)";    # older ,$DIMFLAG)";
		$k = "rct($nrct,$DIMFLAG)";
	}
      } # end of -w=0 bit
	printall("END  Rate $oldrate\n -------------------\n");
	my $retval = $k;   # Returns rate coeffient string

} # end of sub define_rates
########################################################################
sub process_lhs {
# Processes left hand side of equation and creates or adds to
# string $LOSS =  .." for fortran code output, where $LOSS is
# the user-defined loss variable, e.g. L or L(iq)
#######################################################################
	my($i,$other) ;
	my($nlhs) ;
	$nlhs = @lhs;  #Length of lhs array

	my $flux = "$rate_label[$neqn]" ;    # Initialise production flux for process_rhs

	for($i=0;$i<$nlhs;$i++){

#             check for e.g. fgas()*diacid  and that species exists !

		$lhs[$i] = check_multipliers( $lhs[$i] ) ;
		printall("CHECK_MATCH FROM LHS $lhs[$i]\n");
		check_match( $lhs[$i] ) ;

		if ( $loss{$lhs[$i]} ){
			$loss{$lhs[$i]} .=
                          ( $CONTLINE . "   + $rate_label[$neqn]" );
		} else {
			$loss{$lhs[$i]} =
                         "\n      $LOSS =" . $CONTLINE . "     $rate_label[$neqn]" ;
		}
		if ( $nlhs == 2 ) {
			$other = ($i==0) ? 1 : 0 ;
                	$loss{$lhs[$i]} .=  "* $CPNEW($lhs[$other] $CPDIM)" ;
		}
		$flux .=  " * $CPNEW($lhs[$i] $CPDIM)" ;
	}

	if ( $nlhs > 2 ) { die "ERROR - too many left hand side terms:\n",
                                   "  line:  $line\n NLHS= $nlhs \n @lhs " };

	print "fff: $nlhs FFF  $flux\n";
	my $retval = $flux ;

} # end of sub set_losses
#######################################################################
sub check_match {
#                 Checks that the species from lhs or rhs exists in
#                 the species array
#######################################################################
	my( $test_spec) = @_ ;   # @_ is argument to subroutine

	printall("CHECKING MATCH TEST SPEC $test_spec\n");
	# Urghhhhh, but from Prog. Perl, p272
	#my $n = $nspecies[$tot] - 1;
	my $n = $nspecies[$tot]; # J08 fix...
	die "N TOO LOW " if $n < 0;
        my @known =  @{ $species[$tot] }[ 1 .. $n ];

	foreach my $s ( @known ) {
	   return ("species found") if  $test_spec eq $s ;
	}
	# Shouldn't get here...
	my $nn = 0;
        foreach my $s ( @known ) {
	   $nn ++;
	   print "testing:${test_spec}XX against known N$nn Spec:${s}XX\n";
	}
	printall("FAILING FOR $test_spec\n (Could also be missing ;)");
	die "Species $test_spec not found at line $linenum   !!!! \n";
} # end of sub check_match
#######################################################################
sub count_tracer {
##                 Checks that the tracer exists in
##                 the tracer array, if not adds it, and counts atoms
##                 Assumes tracer name and formula are identical. (see docs)
#######################################################################
        my( $test_spec) = @_ ;   # @_ is argument to subroutine

	# Urghhhhh, but from Prog. Perl, p272
	my $n = $nspecies[$tot] - 1;
# if ( $n < 0 ) { print "ZEROCOUNT_TRACER S$test_spec\n"; }
        my @known =  @{ $species[$tot] }[ 0 .. $n ];
        foreach my $t ( @known ) {  #  @species[$tot] ) {
          next unless defined $t;
   	#print "COUNT_TRACER T$t N$n S$test_spec\n";
   	   return ("already known") if ($test_spec eq $t);
        }
#print " NEWTRACER$test_spec\n";
	# If we got to here, we have a new tracer
        count_atoms($test_spec,$test_spec);
        printall ("New tracer, $test_spec,  processed, mol. wt is $molwt{$test_spec} \n");
} # end of sub count_tracer

#######################################################################
sub process_rhs {
## Processes right hand side of equation and creates or adds to
## string "P =  .." for fortran code output. The user-defined variable
## $PROD can be e.g. "P" or "P(iq)" or whatever.
#######################################################################
	my ($spec,$factor,$flux) = @_ ;
	my ($i) ;

	if ( $factor != 1 ){ $flux = $factor . "*" . $flux }

	if ( $prod{$spec} and $factor > 0 )
	{
		$prod{$spec}    .= ( $CONTLINE . "   + $flux" )  ;
	}
	elsif ( $prod{$spec} and $factor < 0 )
	{
		$prod{$spec}    .= ( $CONTLINE . "    $flux" )  ;
	}
	else
	{
		$prod{$spec}    = "\n      $PROD = " . $CONTLINE . "     " .  $flux ;
	}

} # end of sub process_rhs

#######################################################################
sub printsplit {
# Splits up lines after 15 newlines. Important assumption is that
# we can use the first line (e.g. P= &) as a pattern for the new
# continuation lines (givving, e.g. P = P + )
 my($copy) = shift;
 my $n = 0;
 my $newline;
print LOG "PRINTSPLIT start N$n for $copy\n" if $DEBUG;

 if ( $copy =~ /\s*P\s+=/ ) {
	$newline = "      P = P & \n";
 } elsif ( $copy =~ /\s*L\s+=/ ) {
	$newline = "      L = L & \n";
 } else {
 	die "ERROROROOROROORO in printsplit\n";
 }
print LOG "PRINTSPLIT NEWLINE N$n  $newline\n" if $DEBUG;


 while ( $copy =~ /\n/g ) {
	$n++ ;
	 #print LOG "COPYADD $n\n";
 }
 my $lastn = $n;
print LOG "PRINTSPLIT LASTN N$n\n" if $DEBUG;

 $n = 0;
while ( $copy =~ /\G(\n)/g ) {
# while ( $copy =~ /(\n)/g ) {
	print LOG "TEST $copy NPRODLINE N$n MATCH$1ENDMATCH\n";# if $DEBUG;
	my $n15 = $n%15;
	print LOG "N$n N15$n15 LAST$lastn\n";# if $DEBUG;
 	if ( $n >0 && $n%15 == 0  && $n != $lastn  ) {
		$copy =~ s/$1/MARKER/ ;
		print LOG "NPRODLINE N$n for $copy\n";# if $DEBUG;
	}
	$n++;
 }
 $copy =~ s/MARKER/$newline/g ;
 return $copy;
} # sub end printsplit
#######################################################################
#
#######################################################################
sub output_prod_loss {
# Writes out P and L terms to GenOut_reactions.inc, together with any head and tail
# text for the fortran code, e.g. do loops etc, plus the call to 2step.
# Note that if a prod or loss term is empty the following code
# will automatically print it with zero instead.
######################################################################
	my( $eqntext, $noloss, $noprod);

	for( my $i=1; $i<=$nspecies[$tot]; $i++){
        open(OUTFILE,">GenOut_Reactions1.inc") if $i==1;
        open(OUTFILE,">GenOut_Reactions2.inc") if $i==$SlowSpec;
		print OUTFILE "\n!-> $species[$tot][$i] \n";
		print OUTFILE $LOOP_HEADTEXT;

		# Default equan with both prod & loss:

		$eqntext = $CHEMEQN_TEXT ;    # Text to call 2-step
 		$noloss = $noprod = 0 ;
                if ( not defined ($loss{$species[$tot][$i]}) ) {
                    $loss{$species[$tot][$i]} = "      ! $LOSS = 0.0\n";
		    $eqntext = $CHEMEQN_NOLOSS ;  # Text to call 2-step
		    $noloss = 1;
		}
                if ( not defined ($prod{$species[$tot][$i]}) ) {
                    $prod{$species[$tot][$i]} = "      ! $PROD = 0.0\n";
		    $eqntext = $CHEMEQN_NOPROD ;  # Text to call 2-step
		    $noprod = 1;
		}
                if ( $noloss && $noprod ) {
		    $eqntext = $CHEMEQN_NOTERMS;  # Text to call 2-step
		}
		if( $species[$tot][$i] eq "RO2POOL" ){
		    print "RO2 POOL FOUND \n";
                    $prod{$species[$tot][$i]} = "      ! $PROD = 0.0\n";
                    $loss{$species[$tot][$i]} = "      ! $LOSS = 0.0\n";
		    #$eqntext = "      xnew(RO2POOL) = sum( xnew(RO2_POOL) )!!!!\n\n"   ;  # Text to call 2-step
		    $eqntext = "      xnew(RO2POOL) = sum( xnew(RO2_GROUP) )!!!!\n\n"   ;  # Text to call 2-step
			print "RO2 POOL EQN $eqntext \n";
		}

		foreach my $eqn ( $prod{$species[$tot][$i]}, $loss{$species[$tot][$i]} ){
		print OUTFILE "!DSGC eqn $eqn\n============\n" if $DEBUG;
		  if( $eqn ){
			$eqn = printsplit($eqn);   # For long-lines
			print OUTFILE "$eqn \n" };
		}

		$eqntext =~ s/SPEC/$species[$tot][$i]/g ;
		print OUTFILE $eqntext ;
		print OUTFILE $LOOP_TAILTEXT;
	        close(OUTFILE) if $i == $SlowSpec-1;
	}

	close(OUTFILE);
} # end of sub output_prod_loss
#########################################################################
sub carbon_count_etc {
# Checks the number of atoms on the left and right hand-sides of
# the equation. From  @lhs @rhs @lhtracer @rhtracer
  my(%sumleft, %sumright, $spec);

  my $in_balance = 0 ;
  foreach my $atom ( keys(%atoms) ){
	$sumleft{$atom} = $sumright{$atom} = 0;
	foreach $spec ( @lhs, @lhtracer ) {
	    $sumleft{$atom} += $count{$spec}{$atom} ;
	}
	foreach $spec ( @rhs, @rhtracer ) {
            $sumright{$atom} += ( $count{$spec}{$atom} * $rhfactor{$spec} );
        }
	if( $sumleft{$atom} != $sumright{$atom} ) {
		if( $in_balance == 0 ) {printall ("Check balance!...........\n")}
		$in_balance = 1;
		printall("$atom: LEFT $sumleft{$atom} RIGHT $sumright{$atom}\n");
	}
  }
} # end of sub carbon_count_etc

#########################################################################
sub simplify_rate {
print LOG "Simplifying complex rate expression for line $linenum
           rate: $rate \n" ;
	my ($ratestring) = @_ ;
	my ($const,$nconst,$nvar) = ( 1.0,0,0) ;
	my $var = "" ;

	return if ( $ratestring =~ /\+/ ) ; # Too many posibilities for me!
	return if ( $ratestring =~ /\*\*/ ) ; # Too many posibilities for me!
	return if ( $ratestring =~ /^aqrc*/ ) ; # Too many posibilities for me!
	return if ( $ratestring =~ /troe/i ) ; # Too many posibilities for me!

	my ( @ratebits ) = split(/[\*]/,$ratestring);
	foreach my $r ( @ratebits ) {
		if ( is_float($r) ) {
			$const = $const*$r  ;
		} else { # must be variable..
			$var .= "*$r" ;
 		}
		print "RATEBIT is $r => Constants: $const Variables: $var  \n" ;
	}
	$rate = $const . $var ;
	$rate =~ s/^ 1 \* //x;  # Remove initial "1*"
	printall ("SIMPLIFIES  $ratestring  -- to ->  \n") ;
}
#########################################################################
sub check_multipliers {
# checks for e.g. fgas(xxx)*diacid. Only copes with one
# multiplier so far.
	my ($lhsterm) = @_ ;
	my ($mult,$spec);
	return $lhsterm  if ( not $lhsterm =~ /\*/ );

        print LOG "Checking for multipliers for line $linenum :
           rate: $rate", "."x50, "\n";

	( $mult, $spec ) = split(/\*/,$lhsterm) ;
	printall("CHECK_MULTIPLIERS FOR $spec \n");
	check_match($spec);
	$rate_label[$neqn] .= "*$mult" ;
	printall("CHECK_MULTIPLIERS FOUND $mult \n");
	return $spec ;    # returns just the species name now.

} # end of check_multipliers
#########################################################################
## =====================================================================
# Writes to both STOUT and LOG file...
sub printall {
        my($string) = @_ ;
	print "$string";
        print LOG "$string";
}

# =====================================================================
# Number checks : is_float, is_integer
# For an amusing and useful discussion of these, see Perl's FMTEYEWTK
# (far more than even you ever wanted to know) page at:
# ftp://ftp.ccs.neu.edu/net/mirrors/ftp.funet.fi/pub/languages/
#       perl/CPAN/doc/FMTEYEWTK/is_numeric.html
# or, maybe faster, through http://www.perl.com
#
# Subroutine to determine if input is a float.
sub is_float {
    local $^W = 0;     # switches off some unwanted messages
    $_[0] =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;
}
# Subroutine to determine if input is an (positive) integer.
sub is_integer {
    local $^W = 0;     # switches off some unwanted messages
    $_[0] =~ /^\d+$/;
}
# =====================================================================
## Subroutine for max
#sub max {
#my($a,$b) = @_;
#return ($a>=$b) ? $a : $b ;
#}
########################################################################
#
#
#
## Output routines: Now for F
#########################################################################
sub print_species {
#
#............... print_species ..........................
#   outputs fortran-format files GenOut_adv_ml.inc
#   which look like, e.g.:
#
#    Genspec_adv_ml.inc
#        NSPEC_ADV=66
#        integer
#    &   OP, OH, O3, HO2
#        parameter( OP=1)      ! Weight = 16
########################################################################
  my ( $i, $top, $aero_txt );
  my ( %txt )       = ( $tot => "tot",  $adv => "adv", $shl => "shl" );
  my ( %nspec_txt ) = ( $tot => "NSPEC_TOT",
                        $adv => "NSPEC_ADV",
                        $shl => "NSPEC_SHL" );
  my($cnum,$nnum,$snum);  # Carbon, nitrogen and sulphur numbers.
  my($module, $end_of_line, $Use);

  my %ixlab     = ( $tot => "",
                    $adv => "IXADV_",
                    $shl => "IXSHL_"    ) ;

  foreach my $s ( $tot, $adv, $shl ) {
    next if $nspecies[$s] == 0 ;
    print "PROCESS $s TXT $txt{$s} NSPEC $nspecies[$s] \n";

    $module = "ChemSpecs_$txt{$s}_ml";
    open(F,">GenOut_$module.inc");
    $Use = "none";
    start_module($module,\*F,$Use);
    print F "!+ Defines indices and NSPEC for $txt{$s} : $description{$s}\n\n";

    my $nspectxt = "module $module\n
           integer, public, parameter :: $nspec_txt{$s}=$nspecies[$s] \n";

    if ($s == $tot && $naerosol == 0) {   # aerosols kept in totals file
      $first_semivol=-999;
      $last_semivol=-999;
    }
    if ($s == $tot) {   # aerosols kept in totals file
      $aero_txt    = "\n  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=$naerosol,   &!   Number of aerosol species
                FIRST_SEMIVOL=$first_semivol, &!   First aerosol species
                LAST_SEMIVOL=$last_semivol     !   Last  aerosol species  \n\n";
    } else {
      $aero_txt    = "\n";
    }

    # 1 ) nspec and species numbers
    print F "!   ( Output from GenChem, sub print_species ) \n";
    print F "\n  $INTPUB $nspec_txt{$s} = $nspecies[$s] \n $aero_txt";

    for($i=1;$i<=$nspecies[$s]; $i++){
      if ( $i == 1 || $i%10 == 0 ) { # Declare integer
        print F	"\n\n  $INTPUB  & \n" ;
        $end_of_line = "    " ;
      }
      printf F  "$end_of_line";
      my $n = $speciesmap[$s][$i] ;   # Gets n for totals list
      printf F " $ixlab{$s}%-12s=%4d", $species[$s][$n], $i ;
      print " S $s IXLABEL $ixlab{$s} SPEC $species[$s][$n] I $i N $n\n";
      $end_of_line = "   &\n  , ";
    }
    print F "\n\n $HLINE  end module $module\n";
    close(F);
  } # loop over $s


  $module = "ChemChemicals_ml";
  $Use    ="
use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
use ChemSpecs_shl_ml, only: NSPEC_SHL
use ChemSpecs_adv_ml, only: NSPEC_ADV";
  open(F,">GenOut_$module.inc");
  start_module($module,\*F,$Use);

  print F <<"END_CHEMSTART";
!/--   Characteristics of species:
!/--   Number, name, molwt, carbon num, nmhc (1) or not(0)

public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

type, public :: Chemical
     character(len=20) :: name
     real              :: molwt
     integer           :: nmhc      ! nmhc (1) or not(0)
     integer           :: carbons   ! Carbon-number
     real              :: nitrogens ! Nitrogen-number
     integer           :: sulphurs  ! Sulphur-number
     real              :: CiStar    ! VBS param
     real              :: DeltaH    ! VBS param
endtype Chemical
type(Chemical), public, dimension(NSPEC_TOT), target :: species
type(Chemical), public, dimension(:), pointer :: &
  species_shl=>null(),&             ! => species(..short lived..)
  species_adv=>null()               ! => species(..advected..)

contains
subroutine define_chemicals()
!+
! Pointers to short lived and advected portions of species
!
  species_shl=>species(1:NSPEC_SHL)
  species_adv=>species(NSPEC_SHL+1:NSPEC_SHL+NSPEC_ADV)
!+
! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
! array, using indices from total list of species (advected + short-lived).
!                                           MW  NM   C    N   S  C*  dH
END_CHEMSTART

  for($i=1;$i<=$nspecies[$tot]; $i++){
    my $spec = $species[$tot][$i];
    $cnum = $count{$spec}{"C"};   # to save interpolating inside string!
    $nnum = $count{$spec}{"N"};   # to save interpolating inside string!
    $snum = $count{$spec}{"S"};   # to save interpolating inside string!
    printf F
    "    species(%-12s) = Chemical(\"%-12s\",%9.4f,%3d,%3d,%4d,%3d,%8.4f,%7.1f ) \n",
      $species[$tot][$i], $species[$tot][$i],$molwt{$spec}, $nmhc{$spec},
      $cnum, $nnum, $snum, $CiStar[$i], $DeltaH[$i];
    print "SPECF ",
      $species[$tot][$i], $species[$tot][$i],$molwt{$spec}, $nmhc{$spec},
      $cnum, $nnum, $snum, $CiStar[$i], $DeltaH[$i], "\n";
  }
  print F "  end subroutine define_chemicals\n";
  print F "end module $module\n $HLINE";
  close(F);
} # end of sub write_species

#########################################################################
#sub print_speciesmap {
#
#     my $module = "ChemSpecs_maps_ml";
#     my $Use  = "
#       use ChemSpecs_adv_ml, only : NSPEC_ADV
#       use ChemSpecs_shl_ml, only : NSPEC_SHL
#     ";
#
#     open(F,">$module.inc");
#     start_module( $module, \*F, $Use);
#
##	print F "
##  !/--   Mapping indices: from advected species to total index
##  !/--   or short-lived (shl) to total.
##  !/ ...... ..   ( from GenChem )
##        \n";
#
#	foreach my $p ( $adv, $shl ) {
#		printf F "   integer, parameter, public, &\n       dimension($nspec_txt{$p}) :: $description{$p} = (/";
#          	for my $i ( 1..$nspecies[$p] ) {
#			my $comma = "," ;
#			if ( $i == 1  ) {
#				$comma = " &\n      " ;
#			} elsif ( $i%10 == 0 ) {
#				$comma = " &\n      ," ;
#			}
#	 		printf F "%1s %3d", $comma, $speciesmap[$p][$i] ;
#		}
#		printf F "  /)\n\n";
#	}
#        print F " end module $module\n $HLINE";
#	close(F);
#} # end of sub
##################################################################################
sub start_module {   # Simply writes start of  module with horizontal lines and
  my($module) = shift ;
  my($F)       = shift ;
  my($UsedVariables)  = shift ;
  my $EndBit ;

  if($UsedVariables =~ "none") {
    $UsedVariables = "";
    $EndBit       = "implicit none";    # Not needed if no use statements
  } else {
    $EndBit       = "implicit none\nprivate";    # Not needed if no use statements
  }

   print $F <<"MODSTART";
$TLINE
module $module
$HLINE
$UsedVariables
$EndBit
MODSTART
#!/ ...... ..   ( from GenChem )
} # end of sub start_module
########################################################################
sub print_rates {
	my ($rctype,$nrc,$p_rc,$p_rctext) = @_ ; # assign subroutine arguments
       #e.g.  rct,  44,  array, array-text
        my(@rc)     = @$p_rc ;
        my(@rctext) = @$p_rctext ;
	my ($i) = 0;
	my ($tmpout);       # used for storing rate string
#######################################################################

  # Now we print out the rate coefficients. Note that we check for lines > 120
  # characters. If found, we split these at each "*" (crude, but effective
  #  so far!)

  # We also follow the unimod convention that the real rate-coeff.
  # is given by say rcvt(1) = 1.0e-34*exp(34.0/T), whereas the
  # rate used in the equations is denoted rct(1), probably obtained
  # by a tabulation in tablerate. We thus swap text rc to rcv in the
  # the definition and tmpout field.

  my $defrc = $rctype ;
  my $Nrctype = "N" . uc($rctype);
  my $module = "ChemRates_" . $rctype . "_ml" ;
  open(F,">GenOut_$module.inc");
  start_module( $module, \*F, $UsedVariables{$rctype} );

  print F "
  !+ Tabulates $description{$rctype}

    public :: set_${rctype}_rates

    integer, parameter, public :: $Nrctype = $nrc   !! No. coefficients

    real, allocatable, save, public, dimension(:,:) :: $defrc

" ;

  if ( $nrc > 0 ) {
        print F "  contains\n  !------------------------------------\n";
        print F "  subroutine set_${rctype}_rates() \n";
        print F  "     logical,save::first_call=.true.\n";
         
	#A08 print F "  subroutine set_${rctype}_rates($params{$rctype}) \n";
	#A08 print F "  $arguments{$rctype}\n";
	#ESX if ( $rctype =~ /misc/ && $nrctroe > 0 ) {
	if ( $nrctroe > 0 ) {
		print F  "     real, dimension(KCHEMTOP:KMAX_MID) :: log300divt, logtdiv300\n";
		print F  "       log300divt(:) = log(300.0*tinv(:))\n";
		print F  "       logtdiv300(:) = log(temp(:)/300.0)\n";

	}
        print F  "     if(first_call)then\n";
        print F  "       allocate($defrc($Nrctype,$Krange))\n";
        print F  "       $defrc=0.0\n";
        print F  "     endif\n";
                                      

        foreach ( @rctext ) {
		if($_){
			$tmpout = "       $_ = $rc[$i] \n" ;
			my $split_at; # = "\*" ;
			if ( length( $tmpout ) > 80 ) {
			        $split_at = "[\,]" ;
			    # Need to avoid split at comma in e.g. rcmisc(12,:)
				$tmpout =~ s/($split_at\w)/$CONTLINE    $1/g ;
			} # And try again if needed...
			#if ( length( $tmpout ) > 80 ) {
			#        $split_at = "[\*]" ;
			#	$tmpout =~ s/($split_at\w)/$CONTLINE    $1/g ;
			#}
			#$tmpout =~ s/rc/rc_/g ;    # Change rct to rc_t
			$tmpout =~ s/\,iq\)/\)/ ;    # Remove dimension flag here
#			$tmpout =~ s/(XT|xt)/t/g ;    # Change XT  to t
			$tmpout =~ s/(XT|xt)/temp/g ;    # Change XT  to t
			print F  $tmpout ;
	  	}
		else {
			print "IS SS@_ THIS NEEDED?!!!!\n";
		}
	  	$i++ ;
        }
	     print F  "       first_call=.false.";
	     print F  "\n  end subroutine set_${rctype}_rates\n";
  }
  print F  "end module  $module\n";
  close(F);

} # end of sub print_rate
###############################################################################
 sub process_alldep {
  my $ddep = shift;
  my $wdep = shift;
  my $adv  = shift;
  my $found = 0;
  #ok: print "CHECKING ALLDEP DDEP$ddep WDEP$wdep ADV$adv : @species_tot\n";

  $ddep = expand_shorthands($ddep);
    printall( "CDEP_EXPANDED ADV:$adv CDryDEP:$ddep\n");
  $wdep = expand_shorthands($wdep);
    printall( "CDEP_EXPANDED ADV:$adv CWetDEP:$wdep\n");

  # 2010 use group hashes

  unless ($ddep eq "-") {
    my $comma = ($nddep<1)?" ":",";
    $ddep =~ s/DRY_//;
    $ddep_txt .= "      $comma depmap( IXADV_$adv, CDDEP_$ddep, -1) & \n";
    $nddep += 1;
  }
  unless ($wdep eq "-") {
    my $comma = ($nwdep<1)?" ":",";
    $wdep =~ s/WET_//;
    $wdep_txt .= "      $comma depmap( IXADV_$adv, CWDEP_$wdep, -1) & \n";
    $nwdep += 1;
  }
  printall("CDEP $ddep $wdep ADV $adv NDDEP $nddep  NWDEP $nwdep LINE".
    (join '',grep{defined;} @_)." \n");
}

###############################################################################
sub printmap_dep {
  printall( "ENTERING CDEP print $nddep NWDEP $nwdep \n");
  open(DDEP,">GenOut_DryDep.inc") or die "FAIL AllDEP\n"; # if $nddep >  0;

  print DDEP "  integer, public, parameter ::  NDRYDEP_ADV  = $nddep\n";
  print DDEP "  type(depmap), public, dimension(NDRYDEP_ADV), parameter:: DDepMap= (/ &\n";
  print DDEP "$ddep_txt";
  print DDEP "  /)\n";
  close(DDEP);

  open(WDEP,">GenOut_WetDep.inc") or die "FAIL AllDEP\n"; # if $nddep >  0;
  print WDEP "  integer, public, parameter ::  NWETDEP_ADV  = $nwdep\n";
  print WDEP "  type(depmap), public, dimension(NWETDEP_ADV), parameter:: WDepMap= (/ &\n";
  print WDEP "$wdep_txt";
  print WDEP "  /)\n";
  close(WDEP);
}

sub printmap_ext {
  my $nkey=keys %extinc;
  my $outline="";
  printall( "ENTERING CEXT print $nkey \n");
  foreach my $k (sort keys %extinc) {
    my ($cext,$mode)=split(/[:;\-\/|]/,$extinc{$k});
        $mode="Wet" unless $mode;               # WET_MODE by default
    $outline.=($outline?",":"")."&\n"
            .sprintf "  ExtEffMap(%-16s,CEXT_%-4s,%3s_MODE)",$k,$cext,uc($mode);
  }
  open(CEXT,">GenOut_AerExt.inc") or die "FAIL AerExt\n";
  print CEXT 
    "integer, public, parameter :: NUM_EXT = $nkey\n"
   ."type(ExtEffMap), public, dimension(NUM_EXT), parameter :: ExtMap=(/$outline/)\n";
  close(CEXT);
  # define AOD group, for some level of backguard compatibility
  @{$grp{'AOD'}}=(sort keys %extinc);  
}
###############################################################################
sub print_groups {
  my $module = "ChemGroups_ml";
  open(GROUPS,">GenOut_$module.f90");
  my $ngroups  = 0;
  my $groupsub = "\n${HLINE}contains\nsubroutine Init_ChemGroups()\n\n";
  my $Use = "use ChemSpecs_tot_ml  ! => species indices\n"
           ."use OwnDataTypes_ml   ! => typ_sp";
  # implement later
  #   use OwnDataTypes_ml, only : gtype  ! => for group defs";
  start_module($module,\*GROUPS,$Use);
  print GROUPS "! Assignment of groups from GenIn.species:\n";
  print GROUPS "public :: Init_ChemGroups\n\n";

  print "RRRRRRRRRRRRRRRRRRRRRRRRRR ; \n";

  # groups assigned from GenIn.species
  # Check first for empty, e.g. SS can be empty
  foreach my $g (keys %grp) {
    my $N = @{$grp{$g}};
    print "DELETE KEY $g\n" if($N<1);
    delete ($grp{$g})       if($N<1);
  }
  foreach my $g (keys %grp) {
    my $N = @{$grp{$g}};
    my $outline = join(",",@{$grp{$g}});
    print "  TESTSS $g  N$N group is: @{$grp{$g}}\n";
    print GROUPS "integer, public, target, save, dimension($N) :: &\n"
                ."  ${g}_GROUP = (/ $outline /)\n\n";
    $ngroups ++;
    $groupsub .="  chemgroups($ngroups)%name=\"$g\"\n"
               ."  chemgroups($ngroups)%ptr=>${g}_GROUP\n\n";
  }

  print GROUPS "\n! ------- RO2 Pool     species ------------------\n";
  my $N = @ro2pool;
  my $outline = join(",",@ro2pool);
  if($N==0) {   # No RO2 pools defined with older style
    $outline = -99 ; # if  $old_style  ;
    $N = 1;
  }
  print GROUPS "integer, public, parameter :: SIZE_RO2_POOL      = $N\n";
  print GROUPS "integer, public, parameter, dimension($N) :: &\n"
              ."  RO2_POOL      = (/ $outline /)\n";

  print GROUPS "type(typ_sp), dimension($ngroups), public, save :: chemgroups\n\n";
  print GROUPS $groupsub;
  print GROUPS "endsubroutine Init_ChemGroups\n $HLINE";
  print GROUPS "endmodule $module\n $HLINE";
  close(GROUPS);
}
###############################################################################
 sub print_femis {
        my $Nemis = @emis_files;
	printall( "ENTERING FEMIS print $Nemis\n");
        open(EMIS,">femis.defaults") or die "FAIL FEMIS\n";
	print EMIS "Name  $Nemis  ";
	foreach my $e ( @emis_files ){
		printf EMIS "%10s", lc($e);
	}
	print EMIS "\n 28    0 ";
	foreach my $e ( @emis_files ){
		printf EMIS "%10.1f", 1.0 ;
	}
	print EMIS "\n ";
        close(EMIS);
}
###############################################################################
 sub print_emislist {
	my ( @emis) = @_;
	printall( "ENTERING EMISLIST \n");
        open(EMIS,">CM_emislist.csv") or die "FAIL FEMIS\n";
	print EMIS  join(",",@emis );
        close(EMIS);
}
###############################################################################
 sub print_emisstuff {
	my ( $nam, @emis) = @_;
        my $MaxLen; # Was 12
	print "EMISstuff NAM $nam FILES @emis\n";
        my $Nemis = @emis;
        $MaxLen = 0; # Was 12
	foreach my $e ( @emis ){
		my $length=length($e);
		$MaxLen = $length if $length > $MaxLen ;
	}
	printall( "ENTERING EMIS print $Nemis\n");
        open(EMIS,">GenOut_Emis$nam.inc") or die "FAIL EMIS_$nam\n";

	print EMIS "  integer, parameter, public ::  NEMIS_$nam  = $Nemis\n";
	print EMIS "  character(len=$MaxLen), save, dimension(NEMIS_$nam), public:: &\n";
	print EMIS "      EMIS_$nam =  (/ &\n";
	my $comma = "";
	foreach my $e ( @emis ){
		printf EMIS "%12s \"%-${MaxLen}s\" &\n", $comma,$e;  # uc($e);
		$comma = ",";
	}
	print EMIS " /)\n ";
        close(EMIS);
}
###############################################################################
sub print_soa {
	my ( $typ, @soa ) = @_ ; # assign subroutine arguments
	my $N = @soa;
	my $out = join( ', ', @soa );
	return unless $N > 0;
	print SOA "\n!$typ : \n";
	print SOA "   integer, public, parameter :: NUM_$typ = $N\n";
	print SOA "   integer, dimension($N), public, parameter :: &\n";
	print SOA "     $typ = (/  $out /)\n" ;
}
###############################################################################
sub expand_shorthands {
	my $s  = shift;
	foreach my $shorthand (keys(%shorthand)){
		if ($s =~ /\b $shorthand \b/x ) {
			$s =~ s/$shorthand/$shorthand{$shorthand}/g;
		}
	}
	return $s;
} # end of expand_shorthands(
#########################################################################
sub process_emisfile {
	my $arg = shift;
	my ( $emis, $files, $rate );
	if( $arg =~ /:/) { # e.g. emisfiles=sox,nox,nh3
		( $emis, $files ) = split(/:/,$arg);
		my @files = split(/,/,$files);
		my $found = 0;
		foreach my $f ( @files ) {
		   foreach my $e ( @emis_files ) {
			$found = 1 if $f eq $e ;
		   }
		   push( @emis_files, lc($f) ) unless $found;
		}
		$rate =  1; #"emisxxxx($f,k)";
		#printall( "RCEMIS RATE: $rate \n" );
		printall( "EMISFILES EMISF @emis_files\n" );
	}
	return $rate;
}
#########################################################################
sub process_emis {
	my $arg = shift;
	my ( $rcemis, $spec, $rate );
	if( $arg =~ /:/) {
		( $rcemis, $spec ) = split(/:/,$arg);
		my $found = 0;
		foreach my $e ( @emis_specs ) {
			$found = 1 if $spec eq $e ;
			printall( "RCEMIS FOUND: $spec \n" ) if $found;
		}
		push( @emis_specs, $spec ) unless $found;
	        $rate = "rcemis($spec,k)" unless $found ; # Only add if not already set

print "PROCESS_EMIS $arg  $rcemis : $spec \n";
		if ( $rcemis eq "RCBIO" ) {
	  	     $rate = "0 !Skip bio rate since rcemis exists" if $found ; # Only add if not already set
		     my $found = 0;
		     foreach my $e ( @bionat_specs ) {
			$found = 1 if $spec eq $e ;
			printall( "RCBIONAT FOUND: $spec \n" ) if $found;
		     }
		     printall( "NEW RCBIONAT FOUND: $spec \n" ) if $found;
		     push( @bionat_specs, $spec ) unless $found;
		}
		printall( "RCEMIS RATE: $rate \n" );
		printall( "RCEMIS EMISF @emis_specs\n" );
		printall( "RCEMIS BIONA @bionat_specs\n" );
	} else {
		$rate = $arg;
		die "No emission file specified for $arg\n";
	}
	return $rate;
}
#########################################################################
sub process_groups {
  my ( $spec, $groups , $dry, $wet )  =  @_ ;
  $groups = uc($groups);
  my @groups = split(/;/,$groups);
  print "TESTGROUP $groups\n";
  foreach my $g ( @groups ) {
    print "TESTG $g\n";
    #$g = process_group_params($spec,$g) if $g =~ /:/;
    #push(@grp[$g],$spec);
    push @{ $grp{$g}},uc($spec);
    foreach my $gg ( keys %grp ) {
      print "  TESTGG $gg group is: @{ $grp{$gg} }\n" if $g eq "SOX";
    }
    if( $wet ne "-" ) { # Creates groups such as WET_OXN, WET_xxxx ... many!
      my $gwet = "WDEP_$g";
      push @{ $grp{$gwet}},uc($spec);
      foreach my $gg ( keys %grp ) {
       # print "  TESTWETGG $gwet $gg group is: @{ $grp{$gg} }\n"; # if $g eq "SOX";
      }
    } # wet
    if( $dry ne "-" ) { # Creates groups such as DRY_OXN, DRY_xxxx ... many!
      my $gdry = "DDEP_$g";
      push @{ $grp{$gdry}},uc($spec);
      foreach my $gg ( keys %grp ) {
       # print "  TESTdry $gdry $gg group is: @{ $grp{$gg} }\n"; # if $g eq "SOX";
      }
    } # dry
  }
} # end of process_groups
#########################################################################
#sub process_group_params {
#	my ( $spec, $in )  =  @_ ;
#   	my @params = split(/:/,$in);
#	  print "TESTG-PARAMS $in\n";
#	  return $params[0];
#	  #push @{ $grp{$g}},uc($spec);
#
#} # end of process_group_params
#########################################################################
#$groups   = $r[ find_rec("Groups",@r ) ];
sub find_rec {
	my ( $txt, @a )  =  @_ ;
	my $n = 0;
	#my $len = @a;
	for($n=0; $n<= @a; $n++){
		last if $a[$n] eq $txt;
	}
	return $n;
}
#########################################################################
sub specs_equal {   # See if names match, ignoring case
    my ( $n1, $n2 ) = @_ ;
    return ( uc($n1) eq uc($n2) );
}
