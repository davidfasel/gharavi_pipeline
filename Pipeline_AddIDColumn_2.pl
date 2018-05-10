#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);
use List::Util qw(first);
use Getopt::Long;  #process command line arguments

my $IDS_HEADER = "Sample";
my $INSERTED_HEADER = "LabID";

my %samples;

my ($MAPPING_FILE, $SAMPLE_IDS_FILE);
my $bSubmitToSeaSeq = 0;

GetOptions (
    'file|f=s'   => \$SAMPLE_IDS_FILE,
    'map|m=s'    => \$MAPPING_FILE,
    'header|h=s' => \$IDS_HEADER,
    'insert|i=s' => \$INSERTED_HEADER
)  or &usage("");

## make a hash of the original sample ID and the corresponding sample ID
open(FILE, $MAPPING_FILE) or die "Unable to open mapping file";

while (<FILE>) {
  chomp;
  my @a = split(/\t/);
  die "Duplicate sample found in mapping file" if exists $samples{$a[0]};
  $samples{$a[0]} = $a[1];
}
close FILE;



## process the annotation file by inserting an additional column with the corresponding ID
open(FILE, $SAMPLE_IDS_FILE) or die "Unable to open main file.";

# print header with new sample ID column inserted
my $header_line = <FILE>;
my @fields = split(/\t/, $header_line);
#my $id_col = first_index{/$IDS_HEADER/} @fields;
my $id_col = first { $fields[$_] =~ /$IDS_HEADER/ } 0..$#fields;
if (not defined $id_col) {
    print STDERR "Error: Couldn't find column named '$IDS_HEADER' to match with mapping file.\n\n";
    exit;
} 

splice(@fields, $id_col, 0, $INSERTED_HEADER);
print join("\t", @fields);

# output the rest of the file inserting the new sample IDs column
while(<FILE>) {
  my @fields = split(/\t/);
  my $s_ids = $fields[$id_col];
  my @ids = split(/,/, $s_ids);
  #@ids = map {s/\(.+//; $_} @ids;
 
  my @new_ids;  
  for my $id (@ids) {
    $id =~ s/(\(.+)//;
    my $geno = $1 || "";
    
    my $new_id = $samples{$id} || ".";
    push(@new_ids, $new_id . $geno);
  }
  
  my $new_ids_string = join(",", @new_ids). ",";
  
  splice(@fields, $id_col, 0, $new_ids_string);
  print join("\t", @fields);

}




## usage: 
# Pipeline_AddIDColumn.pl -m MappingList.tsv -f File.tsv -h Sample_ID -i New_ID \
#    > AnnotationFile_withNewIDs.tsv
# Input of MappingList.tsv should be the ID in the annotation file
# followed by the new ID to be added to the annotation file
# -h specifies the column name to be matched  (defaults to 
# -i specifies what to name the new column with the new IDs (defaults to LabID)
