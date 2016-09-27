#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);
use List::Util qw(first);

my $IDS_HEADER = "Sample_ID";
my $INSERTED_HEADER = "LabID";

my %samples;
my $SAMPLE_IDS_FILE = $ARGV[0];
my $ANNOTATION_FILE = $ARGV[1];

## make a hash of the original sample ID and the corresponding sample ID
open(FILE, $SAMPLE_IDS_FILE) or die "Unable to open sample IDs file";

while (<FILE>) {
  chomp;
  my @a = split(/\t/);
  die "duplicate sample found in samples file" if exists $samples{$a[0]};
  $samples{$a[0]} = $a[1];
}
close FILE;



## process the annotation file by inserting an additional column with the corresponding ID
open(FILE, $ANNOTATION_FILE) or die "Unable to open annotation file.";

# print header with new sample ID column inserted
my $header_line = <FILE>;
my @fields = split(/\t/, $header_line);
#my $id_col = first_index{/$IDS_HEADER/} @fields;
my $id_col = first { $fields[$_] =~ /$IDS_HEADER/ } 0..$#fields;
splice(@fields, $id_col, 0, $INSERTED_HEADER);
print join("\t", @fields);

# output the rest of the file inserting the new sample IDs column
while(<FILE>) {
  my @fields = split(/\t/);
  my $ids = $fields[$id_col];
  
  for my $key (keys %samples) {
    #todo: be more precise about an exact match on sample ID, in case some sample IDs are contained within
    # longer IDs (ex. Samp1 and Samp10).  So use word boundaries to create an exact match
    $ids =~ s/\b$key\b/$samples{$key}/i;
  }
  
  splice(@fields, $id_col, 0, $ids);
  print join("\t", @fields);
}




## usage: 
# Pipeline_AddIDColumn.pl SampleList.tsv AnnotationFile.tsv > AnnotationFile_withNewIDs.tsv
# Input of SampleList.tsv should be the sample id in the annotation file
# followed by the sample ID to be added to the annotation file
