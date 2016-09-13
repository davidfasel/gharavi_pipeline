#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Pod::Usage;
use List::MoreUtils qw(first_index);
use Scalar::Util qw(looks_like_number);

my ($INFO_COL, $SAMPLE_COUNT_HEADER, $SAMPLEID_HEADER) = (7, 'Samples_Having_Variant', 'Sample_ID(GT)');
my $CADD_FILTER = 10;
my $POLYPHEN_FILTER = 0.5;
#todo: expose in arguments
my $Include_Only_Heterozygous = 1;      
#my $All_Samples_Have_Minor_Allele = 0; 


if (not @ARGV or grep {/-h$|-help$/} @ARGV) { 
  pod2usage(-message=>"", -verbose=>2, -exitval=>2);
}
open(FILE, $ARGV[0]) or die "File $ARGV[0] not found.";



my $line = <FILE>;
print $line;

# get the index of Genotypes and Samples with Variant columns
my @fields = split "\t", $line;
my $Genotypes_Col = first_index {$_ eq $SAMPLEID_HEADER} @fields;
$Genotypes_Col > -1 or die "$SAMPLEID_HEADER not found in header";

# look at every value in the Samples with Variant column to determine the number of 
# individual in the file
# my $Sample_Count_Col = first_index {$_ eq $SAMPLE_COUNT_HEADER} @fields;
# $Sample_Count_Col++;
# my @sample_counts = split("\n", `cut -f $Sample_Count_Col $ARGV[0] | tail -n +2 | sort | uniq`);
# my $Samples_In_File = max(@sample_counts); 

while (<FILE>) {
    my @fields = split "\t";

    # filter by Heterozygous genotypes
    my @genos = split(",", $fields[$Genotypes_Col]);
    if ($Include_Only_Heterozygous) {
        my $het_count = grep {/0\/1/ or /1\/0/} @genos;
        next if ($het_count != scalar @genos);
    }
    
#     if ($All_Samples_Have_Minor_Allele) {
#         next if ($Samples_In_File != scalar @genos);
#     }
    
    #convert Info field into HASH table
    my %info_hash;
    my $info = $fields[$INFO_COL];
    for (split(/;/, $info)) {
        my ($key, $value) = split("=", $_);
        $info_hash{$key} = $value;
    }
    
    # get the CADD and Polyphen scores
    my ($cadd, $polyphen);
    if (looks_like_number($info_hash{'CADD'})) { $cadd = $info_hash{'CADD'} } 
    if (looks_like_number($info_hash{'PH'})) { $polyphen = $info_hash{'PH'} } 
             
    # set Functional / Deleterious Flags
    my ($isLOF, $isMissense);
    if ($info =~ /splice|[^n]frame.?shift|stop.?gain|stop.?los|start.?gain|start.?los/i) {$isLOF = 1;}
    elsif ($info =~ /missense/i) {$isMissense = 1;}
    
    # Filter by CADD and PolyPhen
    next if ($isLOF && defined $cadd && $cadd < $CADD_FILTER);
    next if ($isMissense && defined $polyphen && $polyphen < $POLYPHEN_FILTER);
      
    print;
}

sub max {
    my $max = shift;
    for (@_) {$max = $_ if $_ > $max;}
    return $max;
}

=head1 SYNOPSIS

Dominant_Analysis.pl input.tsv > output.tsv

=head1 DESCRIPTION

Remove variants if all genotypes are not heterozygous.

This assumes that all Individuals are affected and also approximately the same
coverage for each individual. 

=head1 OPTIONS

 -h, -help   Displays this message

Todo:

 -Relies on unix filters to determine the total number of samples (won't work in Windows). 
 -Doesn't read from STDIN (because file name is needed for unix filters).  
 -May change in the future to allow user to choose which Individuals are affected. 
 -May allow user to choose to filter by Het, all affected, or both (currently both).

=head1 AUTHOR

David Fasel
