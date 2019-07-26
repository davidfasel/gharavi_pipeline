#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(first_index);

my @known_patho_vars = (
  "5-112175211-T-A",  # APC variant from Baylor
  "22-29091207-G-A",  # CHEK2 variant from Baylor
  
);

# convert to hash for easy lookup
my %known_patho_vars = map {$_ => 1} @known_patho_vars;
my $header = <STDIN>;
my @header_fields = split("\t", $header);
my $omim_col    = first_index { /omim/i } @header_fields;
my $pop_max_col = first_index { /MAF_max/i } @header_fields;
my $cadd_col    = first_index { /CADD_ann/i } @header_fields;
my $clin_col    = first_index { /CLINSIG/i } @header_fields;
my $pharma_col  = first_index { /pharma/i } @header_fields;
my $metaSVM_col = first_index { /MetaSVM_pred/i } @header_fields;
my $hgmd_col    = first_index { /HGMD_disease/i } @header_fields;
my $sample_count_col = first_index { /Samples_With_Variant/i } @header_fields;
my $variant_id_col   = first_index { /CHROM-POS-REF-ALT/i } @header_fields;
my $source_col   = first_index { /Emerge$/i } @header_fields;
my $genotype_col = first_index { /GENOTYPES$/i } @header_fields;




if ($omim_col == -1 or $pop_max_col == -1 or $cadd_col == -1 or $clin_col == -1 or $pharma_col == -1) {
    die "couldn't find one of the column headers";
}

chomp $header;
print "$header\tCDSReason\n";

while(my $line = <STDIN>)
{ 
    chomp $line;  
    #exit if $. > 10;
    my @line_fields = split("\t", $line);
    
    my $pop_max = $line_fields[$pop_max_col];
    my $omim = $line_fields[$omim_col];
    my $cadd = $line_fields[$cadd_col];
    my $clinvar = $line_fields[$clin_col];
    my $pharma = $line_fields[$pharma_col];
    my $variant = $line_fields[$variant_id_col];
    my $metasvm = $line_fields[$metaSVM_col];
    my $sample_count = $line_fields[$sample_count_col];
    my $hgmd = $line_fields[$hgmd_col];
    my $source = $line_fields[$source_col];
    my $genotype = $line_fields[$genotype_col];
    
    $cadd    = 0 if $cadd eq ".";
    $pop_max = 0 if $pop_max eq ".";
    
    my $print_line = 0;
    my $reason = '';
    
    # if pharmaco-genomic, output line
    #TODO pharma
#     if ($pharma ne ".") {
#         #$print_line = 1;
#     }
    
    
    if (exists $known_patho_vars{$variant} ) {
      $print_line = 1;
      $reason = 'Known';
    }
    
    # skip if it's common within the cohort
    elsif ($sample_count > 10) {
        next;
    }
    #skip if heterozygous and source says Autosomal Recessive
    elsif (($genotype =~ /^0/ or $genotype =~ /^..0/) and  $source =~ /AR/i ) {
        next;
    }
    
    # Proceed if rare and omim is known and clinvar does not contain 'benign'
    elsif ($pop_max <= 0.001 and $omim ne "." and $clinvar !~ /benign/i){
        $_ = $line;  # just a shortcut so the regex's below apply to the whole line
        
        # now check if HGMD or Clinvar call it pathogenic
        if ($clinvar =~ /patho/i) {
            $print_line = 1;
            $reason = 'Clinvar';
        }
        
        elsif ($hgmd =~ /DM/i and $hgmd !~ /DM\?/) {
            $print_line = 1;
            $reason = 'HGMD';
        }
        
        # if splicing, stopgain/loss, frameshift, inframe-indels, start-gain
        elsif ( /stop.*gain/i || /stop.*los/i || /start_lost/i || 
                /\tframe_?shift/i || /exon_del/i || (/\tsplic/i and $line !~ /region/i) ) {
            $print_line = 1;
            $reason = 'LOF';
        }
        
        # inframe-indels
#         elsif (/codon.*(del|ins)/i || /nonframeshift/i || /cds-indel/i || /coding-near-splice/i) {  
#             $print_line = 1;
#         }
        
        # if very rare and metaSVM predicts damaging (used to be high CADD)
        elsif ($pop_max <= 0.00001 and $metasvm =~ /D/i ) {
            $print_line = 1;
            $reason = 'RareMetaSVM';
        }
    } 
    
    if ($print_line) {
        $line_fields[$source_col] = "CUMC_CDS" if $source eq ".";
        print join("\t", @line_fields) . "\t$reason\n" ;
    }
    
    #last if $. > 100;

}

