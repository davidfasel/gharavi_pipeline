#!/usr/bin/perl
use warnings;
use strict;

my $GENOTYPES_COL = 9;
# my $MRN_COL = 10;
# my $EMERGE_COL = 11;
# my $SAMPLE_COL = 12;
my $SAMPLE_COL = 10;


open(FILE, $ARGV[0]) or die "Can't open $ARGV[0].";
my $header = <FILE>;
print $header;

my @col_headers = split("\t", $header);

     
while (my $line = <FILE>) {
  my @fields = split("\t", $line);
  my @genotypes = split(";", $fields[$GENOTYPES_COL]);
  # my @MRNs = split(",", $fields[$MRN_COL]);
  # my @emerge_ids = split(",", $fields[$EMERGE_COL]);
  my @sample_ids = split(",", $fields[$SAMPLE_COL]);
  
  # for my $value (@MRNs) {$value =~ s/\(...\)//};
  # for my $value (@emerge_ids) {$value =~ s/\(...\)//};
  for my $value (@sample_ids) {$value =~ s/\(...\)//};
  
  for (my $i = 0; $i < @sample_ids; $i++) {
    # print join("\t", @fields[0..8], $genotypes[$i], $MRNs[$i], $emerge_ids[$i], $sample_ids[$i], @fields[13 .. $#fields]);      
    print join("\t", @fields[0..8], $genotypes[$i], $sample_ids[$i], @fields[11 .. $#fields]);

  }

} 