#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);

# Header names, must be updated if header names change;
my $CLINVAR = "CLINSIG";
my $KIDNEY  = "KidDisorder";
my $OMIM    = "OMIM_Disorder";
my $EMERGE  = "Emerge";
my $EMERGE_SNP = "EmergeSNP";
my $MOST_DELETERIOUS = "MostDeleterious";


my $Input_File = $ARGV[0];
my $Known_Genes_File = $Input_File;
my $Deleterious_File = $Input_File;
$Known_Genes_File =~ s/.tsv$/_Kid_Clin_Omim_Emrg.tsv/;
$Deleterious_File =~ s/.tsv$/_Deleterious.tsv/;

open(FILE, $Input_File) or die "Unable to open sample Input file";
open(KNOWN_GENES, "> $Known_Genes_File") or die "Can't open output file $Known_Genes_File";
open(DELETERIOUS, "> $Deleterious_File") or die "Can't open output file $Deleterious_File";



## print the file filtered by Known disorders from Clinvar, etc.
my $header = <FILE>;
my @fields = split(/\t/, $header);

my $clinvar_col = first_index {/$CLINVAR/} @fields;
my $kidney_col  = first_index {/$KIDNEY/} @fields;
my $omim_col    = first_index {/$OMIM/} @fields;
my $emerge_col  = first_index {/$EMERGE/} @fields;
my $emerge_snp_col = first_index {/$EMERGE_SNP/} @fields;
my $most_deleterious_col = first_index {/$MOST_DELETERIOUS/} @fields;

print KNOWN_GENES $header;
print DELETERIOUS $header;

while (my $line = <FILE>) {
    my @fields = split(/\t/, $line);
    if ($fields[$clinvar_col] ne "." || $fields[$kidney_col] ne "." || $fields[$omim_col] ne "." ||
        $fields[$emerge_col] ne "." || $fields[$emerge_snp_col] ne ".")   {
        print KNOWN_GENES $line;
        
    }
    
    # a stricter filter for the Deleterious list from Pipeline_SelectFunctional.pl
    if (not $fields[$most_deleterious_col] =~ /INTRAGENIC|SPLICE_SITE_REGION|ncRNA|non-coding/ ) {
        print DELETERIOUS $line;
    }
}

