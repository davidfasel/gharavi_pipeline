#!/usr/bin/perl
use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(first_index);

# Update these headers if they are changed in Pipeline_StatusPublicDB.pl
my @HEADER_LABELS = qw(Samples_Having_Variant 141_AF  Tot_AF  EUR_AF  AFR_AF  all  afr  eur  MAF_max);

usage() if (@ARGV != 3);
my ($file, $pop, $freq_thresh) = @ARGV;
chomp (my $line_count = `wc -l $file`);

# process the header
open(FILE, "$file") or die "Couldn't open $file";
my $header = <FILE>;
$header =~ /^CHROM/ or die "Invalid header in $file";
my @fields = split("\t", $header);
my $SAMPLES_WITH_VARIANT = first_index {$_ eq $HEADER_LABELS[0]} @fields;
my $DBSNP_FREQ = first_index {$_ eq $HEADER_LABELS[1]} @fields;
my $TOT_ESP    = first_index {$_ eq $HEADER_LABELS[2]} @fields;
my $EUR_ESP    = first_index {$_ eq $HEADER_LABELS[3]} @fields;
my $AFR_ESP    = first_index {$_ eq $HEADER_LABELS[4]} @fields; 
my $TOT_1K     = first_index {$_ eq $HEADER_LABELS[5]} @fields;
my $AFR_1K     = first_index {$_ eq $HEADER_LABELS[6]} @fields; 
my $EUR_1K     = first_index {$_ eq $HEADER_LABELS[7]} @fields; 
print $header;

while (my $line = <FILE>) {
    print STDERR "processing line $. of $line_count \r";
    
    my @fields = split(/\t/,$line);
    my @dbsnp_freq = split(/,/, $fields[$DBSNP_FREQ]);

    # look for the lowest alt allele freq ($i starts at 1, because we skip the first freq which is the ref allele)
    my $dbsnp_altallele_freq = 0;
    for(my $i = 1; $i < @dbsnp_freq; $i++) {
        my $temp = &makeNum($dbsnp_freq[$i]);
        if ($i == 1) {
          $dbsnp_altallele_freq = $temp;
        }
        elsif ($dbsnp_altallele_freq > $temp) {
          $dbsnp_altallele_freq = $temp;
        }
    }

    my $tot_esp_freq = &makeNum($fields[$TOT_ESP]);
    my $eur_esp_freq = &makeNum($fields[$EUR_ESP]);
    my $afr_esp_freq = &makeNum($fields[$AFR_ESP]);
    my $tot_OneK_freq = &makeNum($fields[$TOT_1K]);
    my $afr_OneK_freq = &makeNum($fields[$AFR_1K]);
    my $eur_OneK_freq = &makeNum($fields[$EUR_1K]);
  
    # Don't print the line if no one has the variant
    # TODO: it's possible that the reference variant is actually the rare variant, so we should check
    next if ($fields[$SAMPLES_WITH_VARIANT] == 0);
    
    my $t = $freq_thresh; 
    $pop = uc $pop;
    if    (($pop eq "T") && ($dbsnp_altallele_freq <= $t) && ($tot_OneK_freq <= $t) && ($tot_esp_freq <= $t)) { }
    elsif (($pop eq "E") && ($dbsnp_altallele_freq <= $t) && ($eur_OneK_freq <= $t) && ($eur_esp_freq <= $t)) { }
    elsif (($pop eq "A") && ($dbsnp_altallele_freq <= $t) && ($afr_OneK_freq <= $t) && ($afr_esp_freq <= $t)) { }
    elsif (($pop eq "EA") && ($dbsnp_altallele_freq <= $t) && ($eur_OneK_freq <= $t) && ($afr_OneK_freq <= $t) && ($eur_esp_freq <= $t) && ($afr_esp_freq <= $t)) { }
    else  { next; }
    print "$line";
}
close FILE;
print STDERR "\n";

sub makeNum () {
  my $num = shift;
  return looks_like_number($num) ? $num : 0;
}

sub usage {
    my $msg = "USAGE: ./select_rare_choose.pl query_file Population(T/E/A/EA) Your_MAF_Threshold
               Population: T(otal)\tE(uropean)\tA(frican)\tEA(Europian and African)
               A MAF threshold 1 in 100 is 0.01, 1 in 1000 is 0.001 and so on";
    print STDERR $msg =~ s/^ +//m;
    exit 1;
}