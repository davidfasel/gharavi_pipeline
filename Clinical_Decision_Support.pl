#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(first_index);


my $header = <STDIN>;
my @header_fields = split("\t", $header);

my $omim_col    = first_index { /omim/i } @header_fields;
my $pop_max_col = first_index { /MAF_max/i } @header_fields;
my $cadd_col    = first_index { /CADD_ann/i } @header_fields;
my $clin_col    = first_index { /CLINSIG/i } @header_fields;
my $pharmo_col    = first_index { /pharmo/i } @header_fields;

if ($omim_col == -1 or $pop_max_col == -1 or $cadd_col == -1 or $clin_col == -1 or $pharmo_col == -1) {
    die "couldn't find one of the column headers";
}

print $header;

while(my $line = <STDIN>)
{   
    my @line_fields = split("\t", $line);
    
    my $pop_max = $line_fields[$pop_max_col];
    my $omim = $line_fields[$omim_col];
    my $cadd = $line_fields[$cadd_col];
    my $clinvar = $line_fields[$clin_col];
    my $pharmo = $line_fields[$pharmo_col];
    
    $cadd    = 0 if $cadd eq ".";
    $pop_max = 0 if $pop_max eq ".";
    
    my $print_line = 0;
    
    # if pharmaco-genomic, output line
    #TODO pharmo
    if ($pharmo ne ".") {
        $print_line = 1;
    }
    
    # Proceed if MAF < 0.01 and OMIM != NA
    elsif ($pop_max <= 0.01 and $omim ne "."){
        $_ = $line;  # just a shortcut so the regex's below apply to the whole line
        
        # now check if HGMD or Clinvar call it pathogenic
        #TODO HGMD
        if ($clinvar =~ /patho/) {
            $print_line = 1;
        }
        
        # if splicing, stopgain/loss, frameshift, inframe-indels, start-gain
        elsif ( /stop.*gain/i || /stop.*los/i || /start.*los/i || 
                /frameshift/i || /exon_del/i ) {
            $print_line = 1;
        }
        
        # inframe-indels
        elsif (/codon.*(del|ins)/i || /nonframeshift/i || /cds-indel/i || /coding-near-splice/i) {  
            $print_line = 1;
        }
        
        # if very rare and high CADD 
        elsif ($pop_max <= 0.001 and $cadd >= 15) {
            $print_line = 1;
        }
    } 
    
    print $line if $print_line;

}