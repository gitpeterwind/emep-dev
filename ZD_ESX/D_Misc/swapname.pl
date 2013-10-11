#!/usr/bin/env perl
# A perl script to swap some unwanted variable names to
# something else.
# e.g. below we swap the name nlc to NCANOPY_LAYERS
# and nlcp to NCANOPY_BOUNDARIES
use strict;

my %swaps = (
  'real(dp)' => "real"
);
# done
#    zmid => "zbnd"
#Not done
#   ,dzmid => "dzcntrs"


while (my $line= <>){
  foreach my $old ( keys %swaps ) {
    my $new = $swaps{$old};
    if ( $line =~ / \b $old \b /x ) {
       $line =~s/ \b ($old) \b /$new/gx
    }
  }
  print $line;
}
