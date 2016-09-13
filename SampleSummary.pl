#!/usr/bin/perl
use strict;
use warnings;
use feature 'say';
use IO::Zlib;
use List::MoreUtils qw(first_index);

# output the header
my @HEADER = (
  'NAME',
  'VARIANTS(AFTER QUAL AND COV FILTER)',
  'Rare/Novel',
  'Hets',
  'Homs',
  'Stopgain',
  'Stoploss',
  'Splicing',
  'Nonsynonymous',
  'Frameshift_InDel',
  'Inframe_InDel',
  'Splicing_InDel',
  'LOFs(stopgain+Splice+Frame+SplicingIndel)',
  'Funtional (stop,splice,missense,frameshift and non-frameshift)',
  'Variants in Kidney Genes',
  'Pathogenic in ClinVardb',
  'Private Mutations'
);
say join("\t", @HEADER);

# parse the header of the input file
my $header = <>;
my @head = split "\t", $header;

my $SampleIndex   = first_index{/Sample_ID/} @head;
my $MAFIndex      = first_index{/MAF_max/} @head;
my $KidneyIndex   = first_index{/KidDisorder/} @head;
my $ClinvarIndex  = first_index{/Clinvar/} @head;
#my $Mutation = first_index{} @head;

# do the counts in each line of the input file

# Stop Gain:
#   stop-gained-near-splice
#   STOP_GAINED
#   stop-gained
#   stopgain
#   stop-gained
# 
# Stop Lost:
#   stop-lost
#   STOP_LOST
#   stoploss
#   
# Frameshift Near Splice:
#   frameshift-near-splice
#   
# Splice:
#   splice-acceptor
#   splice-acceptor-variant
#   SPLICE_SITE_ACCEPTOR
#   splice-donor
#   splice-donor-variant
#   SPLICE_SITE_DONOR
#   
# Frameshift:
#   frameshift_deletion
#   frameshift_insertion
#   frameshift-variant
#   frameshift
#   FRAME_SHIFT
# 
# Missense:
#   missense-near-splice
#   missense
#   NON_SYNONYMOUS_CODING
#   nonsynonymous_SNV
#     
# Inframe:
#   CODON_CHANGE_PLUS_CODON_DELETION
#   CODON_CHANGE_PLUS_CODON_INSERTION
#   CODON_DELETION
#   CODON_INSERTION
#   nonframeshift_deletion
#   nonframeshift_insertion
#   nonframeshift_substitution
#   cds-indel

my %SampleHash;
while (<>) {
    print STDERR "Processing line $.\r" if $. % 100 == 0;
    next if not /\S+/;
    my @fields = split(/\t/);
    my @samples = split(/,/, $fields[$SampleIndex]);

    #s/\(.+\)// for @samples;  #remove genotype from sample name
  
    for my $s (@samples) { 
        my $sample = $s;  # create a temp copy of sample before trimming genotype
        $s =~ s/\(.+\)//;  #remove genotype from sample name        
        $SampleHash{$s}{'het'}++ if $sample =~ /\(0\/\d\)|\(1\d\/0\)/;
        $SampleHash{$s}{'hom'}++ if $sample =~ /\(1\/1\)/;
 
        $SampleHash{$s}{'count'}++;
  
        ### see Pipeline_SelectFunctional.pl for list of functional mutations
        if (/stop.?gain/i) {
            $SampleHash{$s}{'stopgain'}++;
        }
        elsif (/stop.?los/i) {
            $SampleHash{$s}{'stoploss'}++;
        }
        elsif (/frameshift-near-splice/i) {
            $SampleHash{$s}{'fssplice'}++;
        }
        elsif (/splice.*accept|splice.*donor/i) {
            $SampleHash{$s}{'splice'}++;
        }
        elsif (/[^n]frame.?shift/i) {
            $SampleHash{$s}{'frameshift'}++;
        }
        elsif (/non.?synonymous|missense/i) {
            $SampleHash{$s}{'missense'}++;
        }
        elsif (/codon_|nonframeshift_del|nonframeshift_ins/i) {
            #removed: |cds-indel
            $SampleHash{$s}{'inframe'}++; 
        }        
  
        #private mutations
        if (@samples == 1  &&  $fields[$MAFIndex] eq '.') {
            $SampleHash{$s}{'private'}++;
        }
        if ($KidneyIndex != -1 && $fields[$KidneyIndex] ne '.') {
            $SampleHash{$s}{'kidney'}++;
        }
        if ($ClinvarIndex != -1 && $fields[$ClinvarIndex] ne '.'  &&  $fields[$ClinvarIndex] =~ /pathogenic/) {
            $SampleHash{$s}{'clinvar'}++;
        }
    }
}


# output the results for each sample to STDOUT
for my $s (sort keys %SampleHash) {
    $SampleHash{$s}{'stoploss'} //= 0;
    $SampleHash{$s}{'inframe'}  //= 0;
    $SampleHash{$s}{'fssplice'} //= 0;
    $SampleHash{$s}{'kidney'}   //= 0;
    $SampleHash{$s}{'frameshift'} //= 0;

    my $lof = $SampleHash{$s}{'stopgain'} + 
              $SampleHash{$s}{'splice'} + 
              $SampleHash{$s}{'fssplice'} +
              $SampleHash{$s}{'frameshift'};

    my $func = $SampleHash{$s}{'stopgain'} + 
               $SampleHash{$s}{'stoploss'} +
               $SampleHash{$s}{'splice'} + 
               $SampleHash{$s}{'fssplice'} +
               $SampleHash{$s}{'frameshift'} + 
               $SampleHash{$s}{'inframe'} +
               $SampleHash{$s}{'missense'};

    say join("\t", 
        $s, "", "",
        $SampleHash{$s}{'het'},
        $SampleHash{$s}{'hom'},
        $SampleHash{$s}{'stopgain'},
        $SampleHash{$s}{'stoploss'},
        $SampleHash{$s}{'splice'},
        $SampleHash{$s}{'missense'},
        $SampleHash{$s}{'frameshift'},
        $SampleHash{$s}{'inframe'},
        $SampleHash{$s}{'fssplice'},
        $lof,
        $func,
        $SampleHash{$s}{'kidney'},
        $SampleHash{$s}{'clinvar'},
        $SampleHash{$s}{'private'},
    )
}

# Stopgain
# Stoploss
# Splicing
# Nonsynonymous
# Frameshift_InDel
# Inframe_InDel
# Splicing_InDel
# LOFs(stopgain+Splice+Frame+SplicingIndel)
# Funtional (stop,splice,missense,frameshift and non-frameshift)
# Variants in Kidney Genes
# Pathogenic in ClinVardb 
# Private Mutations
