#!/usr/bin/env perl
use strict;

my ( $Spec,$adv,$forumula,$MW,$Dry,$Wet,$Extinc,
      $Cstar,$DeltaH,$Empty,$Groups,$Comments );
my  ( $spec,  $adv, $formula, $wt, $c, $comment );

open(SPEC,"<CRI_v2_0old.species");
open(NEW,">CRI_v2_0base.species");
my $Num = 0;
while(<SPEC>){
   $Num++ ;
   print NEW  if $Num < 5;    #/^\*/;
   ( $spec,  $adv, $formula, $MW, $c, $comment ) = split;
   print NEW   "\"#GROUPED_STYLE - fewer columns than Garry's new-style files\"\n" if $Num == 5; # $adv eq "Species";

   next if $Num < 8;
   print NEW "\"\#SLOW\n"  if $spec =~ /.*SLOW.*/;
#print "SPEC $spec ABD $adv\n";
   next if $spec =~ /.*SLOW.*/;
#print "NEW $Num\n";
 $Dry = "xx";
 $Wet = "xx";
 $Extinc = 0;
 $Cstar  = 0;
 $DeltaH = 0;
 $Groups = 0;


 if ( $comment eq "peroxide" ) {
   $comment = "RO2";
   $Dry     = "ROOH";
 } elsif ( $comment eq "peroxy" ) {
   $comment = "RO2";
 } elsif ( $comment eq "aldehyde" ) {
   $Dry     = "RCHO";
 } elsif ( $comment eq "PAN" ) {
   $Dry     = "PAN";
 } elsif ( $comment eq "biogenic" ) {
   $Groups    = "BVOC";
 } elsif ( $comment eq "PM" ) {
   my $pm = "PMf";
   $pm = "PMc" if $spec =~ /PMco/;
   $Dry    = $pm;
   $Wet    = $pm;
 } elsif ( $comment eq "SOx" ) {
   $Groups  = "SOX";
 } elsif ( $comment eq "NOx" ) {
   $Groups  = "OXN";
   $Dry = "NO2"  if $spec eq "NO2";
 } elsif ( $comment eq "NOy" ) {
   $Groups  = "OXN";
   $Dry = "HONO"  if $spec eq "HONO";
   $Wet = "HNO3"  if $spec eq "HONO";
   $Dry = "HNO3"  if $spec eq "HNO3";
   $Wet = "HNO3"  if $spec eq "HNO3";
   $Dry = "PMfN"  if $comment eq "PM";
 }

 $Dry = "PMfS"  if $spec eq "SO4";
 $Wet = "SO4"  if $spec eq "SO4";
 $Dry = "NH3"  if $spec eq "NH3";
 $Wet = "NH3"  if $spec eq "NH3";
 $Groups  = "RDN"  if $spec eq "NH3";
 $Groups  = "RDN"  if $spec eq "NH4_f";
 $Dry = "O3"  if $spec eq "O3";
 $Dry = "SO2"  if $spec eq "H2O2";   # OBS!
 $Wet = "H2O2"  if $spec eq "H2O2";
 $Dry = "HCHO"  if $spec eq "HCHO";   # OBS!
 $Wet = "HCHO"  if $spec eq "HCHO";
 $Dry = "SO2"  if $spec eq "SO2";   # OBS!
 $Wet = "SO2"  if $spec eq "SO2";   # OBS!
 $Dry = "NH3"  if $spec eq "NH3";   # OBS!
 $Wet = "NH3"  if $spec eq "NH3";   # OBS!
 $MW = "xx" if $MW eq "-";
 print NEW "$spec,$adv,$formula,$MW,$Dry,$Wet,$Extinc, $Cstar,$DeltaH,$Empty,$Groups,xx,xx,,$comment!\n";
}

#my ( $Spec,$adv,$forumula,$MW,$DRY,$WET,$Extinc,
#      $Cstar,$DeltaH,$Empty,$Groups,$Comments );
