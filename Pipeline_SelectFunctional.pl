#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(uniq);
use Array::Utils qw(intersect);

my @HEADER = qw(==  MostDeleterious);
my $MISSING = ".";

# summary list from above, ordered from most deleterious
my @Deleterious = qw(
  stop-gained-near-splice
  STOP_GAINED
  stop-gained
  stopgain
  stop-gained
  START_LOST
  
  frameshift-near-splice
  frameshift_deletion
  frameshift_insertion
  frameshift-variant
  frameshift
  FRAME_SHIFT
  EXON_DELETED

  splice-acceptor
  splice-acceptor-variant
  SPLICE_SITE_ACCEPTOR
  splice-donor
  splice-donor-variant
  SPLICE_SITE_DONOR
  SPLICE_SITE_REGION
  coding-near-splice
  coding-unknown-near-splice
  synonymous-near-splice
  non-coding-exon-near-splice
  intron-near-splice
  splicing

  stop-lost	
  STOP_LOST
  stoploss

  CODON_CHANGE_PLUS_CODON_DELETION
  CODON_CHANGE_PLUS_CODON_INSERTION
  CODON_DELETION
  CODON_INSERTION
  nonframeshift_deletion
  nonframeshift_insertion
  nonframeshift_substitution
  cds-indel

  RARE_AMINO_ACID
  NON_SYNONYMOUS_START
  missense-near-splice
  START_GAINED
  missense
  NON_SYNONYMOUS_CODING
  nonsynonymous_SNV
  coding
  codingComplex
  coding-unknown
  
  UTR_5_DELETED

  ncRNA_exonic
  ncRNA_splicing
  non-coding-exon
  nc-transcript-variant
  INTRAGENIC
);  

my @Benign = qw(
  .
  downstream
  DOWNSTREAM
  downstream-gene
  downstream-variant-500B
  exonic
  EXON
  intergenic
  INTERGENIC
  ncRNA_intronic
  ncRNA_UTR5
  intron
  INTRON
  intron-variant
  intronic
  NON_SYNONYMOUS_START
  NON_SYNONYMOUS_STOP
  START_GAINED
  synonymous
  SYNONYMOUS_CODING
  synonymous_SNV
  SYNONYMOUS_START
  SYNONYMOUS_STOP
  synonymous-codon
  upstream
  UPSTREAM
  upstream-gene
  upstream-variant-2KB
  UTR_3_DELETED
  UTR_3_PRIME
  UTR_5_PRIME
  utr-variant-3-prime
  utr-variant-5-prime
  UTR3
  UTR5
  3-prime-UTR
  5-prime-UTR
  unknown
  NONE
);

my %AllMutationTypes;
for(@Deleterious, @Benign) {
    $AllMutationTypes{$_} = 1;
}

open(FILE, $ARGV[0]) or die "Can't open $ARGV[0].";
my $line = <FILE>;
chomp $line;
die "Invalid Header in file $ARGV[0]." if($line !~ /^CHROM/);
      
# get the index of the INFO column
my $InfoIndex = -1;
my @fields = split(/\t/, $line);
for my $header (@fields) {
    $InfoIndex++;
    last if $header eq "INFO";
}

#print the header
say join ("\t", $line, @HEADER);  

while(my $line = <FILE>) {
    chomp($line);
    my @output = ($line);
    my @fields = split(/\t/, $line);
    
    #convert Info field into HASH table
    my %info_hash;
    my @info_field = split(/;/, $fields[$InfoIndex]);
    for my $item (@info_field) {
      my ($key, $value) = split("=", $item);
      $value //= ".";
      $info_hash{$key} = $value;
    }

    my (@mutationSeaSeq, @mutationsAnnovar, @mutationsSNPEFF);
    $info_hash{"Func.refGene"} //= "";
    $info_hash{"ExonicFunc.refGene"} //= "";
    $info_hash{"FG"}           //= "";
    $info_hash{"FD"}           //= "";
    $info_hash{"EFF"}          //= "";
    @mutationsAnnovar = (
      split(/,|\\x3b/, $info_hash{"Func.refGene"}), 
      split(/,|\\x3b/, $info_hash{"ExonicFunc.refGene"})
    );

    push(@mutationSeaSeq, split("[,/]", $info_hash{"FG"})) if $info_hash{"FG"};
    push(@mutationSeaSeq, split("[,/]", $info_hash{"FD"})) if $info_hash{"FD"};

    my @snpEff = split(",", $info_hash{"EFF"});
    for my $transcript (@snpEff) {
        my ($mutationType) = split(/\(/, $transcript);
        push(@mutationsSNPEFF, $mutationType);
    }
    
    my @mutations = (@mutationSeaSeq, @mutationsAnnovar, @mutationsSNPEFF);
    
    #look for an unknown mutation type
    my $unknown = 0;
    for (@mutations) {
      if(not $AllMutationTypes{$_}) {
          say STDERR "Unknown mutation type '$_', so the variant will be kept. Line $.: $fields[$InfoIndex]. ";
          $unknown = 1;
      }
    }
    
    my $most_damaging = &getMostDamaging(\@mutations, \@Deleterious);
    
    # goto next line if we can't find one of the functional mutations we are interested in and none are unknown
    next if not $most_damaging and not $unknown;
  
    push(@output, "==", $most_damaging); 
    say join ("\t", @output)
    
}  # end main while loop
close FILE;


# getMostDamaging (@list, @refList )
# returns the most damaging mutation in list (assuming that refList is ordered by most damaging first)
# If there is an intersection between the list and the reference list,
# then order the list according to the reference list and return the first item.
sub getMostDamaging () {
    my @list = uniq( @{$_[0]} );
    my @inter = intersect(@list, @{$_[1]});
    return $inter[0] if @inter;
    return 0;
}

