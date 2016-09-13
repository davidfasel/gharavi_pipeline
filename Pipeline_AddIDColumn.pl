#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);

my $IDS_COLUMN = 10;  #zero based column indexed 
my $IDS_HEADER = "Sample_ID";
my $INSERTED_HEADER = "LabID";



my %patients;

open(FILE, $ARGV[0]) or die "Unable to open patient IDs file";
while (<FILE>) {
  chomp;
  my @a = split(/\t/);
  die "duplicate seq number found" if exists $patients{$a[0]};
  $patients{$a[0]} = $a[1];
}
close FILE;



open(FILE, $ARGV[1]) or die "Unable to open annotation file.";

#insert column into header
my $header_line = <FILE>;
my @fields = split(/\t/, $header_line);
my $id_col = first_index{/$IDS_HEADER/} @fields;
splice(@fields, $id_col, 0, $INSERTED_HEADER);
print join("\t", @fields);


while(<FILE>) {
  my @fields = split(/\t/);
  my $ids = $fields[$id_col];
  for my $key (keys %patients) {
    $ids =~ s/$key/$patients{$key}/;
  }
  
  splice(@fields, $id_col, 0, $ids);
  print join("\t", @fields);
}


