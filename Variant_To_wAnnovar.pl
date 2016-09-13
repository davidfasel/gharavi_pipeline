#!/usr/bin/perl -lan
use warnings;
use strict;

my $pos_end = $F[1];

# First column must be a digit
if ($F[0] !~ /^\d/) { next}

if (length($F[2]) > 1) {  #deletion
  $F[3] = "-";
  $F[2] = substr($F[2], 1);
  $pos_end = $F[1] + length($F[2]) - 1;
}
elsif (length($F[3]) > 1) {  #insertion
  $F[2] = "-";
  $F[3] = substr($F[3], 1);
}

splice (@F, 2, 0, $pos_end);
print join("\t", @F);
