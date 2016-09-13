#!/usr/bin/perl
use strict;
use warnings;
use feature 'say';
use IO::Zlib;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq first_index);

my $ESP_FILE = "/media/Data/public_databases/ESP6500_V2_Table2.xls";
my ($INFO_COL, $GENOTYPES_HEADER) = (7, 'Sample_ID(GT)');
my $DATA_NA = ".";

my ($File_count) = (0);
my (%genes, %gene_lengths, %Individuals);

# parse the ESP file to get gene lengths 
open(ESP, $ESP_FILE) or die "Couldn't open ESP file.";
while (<ESP>) {
    my @fields = split("\t");
    if (looks_like_number($fields[1])) {
        $gene_lengths{$fields[0]} = $fields[1];
    }
}

# parse the input files
for my $file (@ARGV) {
    my %file_genes;
    
    # determine if the file is gzipped and open it
    my $fh;
    if ($file =~ /\.gz$/) {
        $fh = new IO::Zlib;
        $fh->open($file, "rb") or die "Unable to open: $file";
    }
    else {
        open($fh, '<', $file) or die "Unable to open: $file";
    }
    
    # get the name of the family based on the file name
    (my $fam = $file) =~ s|.+/||; #remove path
    $fam =~ s/\..+$//;  #remove extensions
    
    # get index of Genotypes col
    my $header = <$fh>;
    my $Genotypes_Col = first_index {$_ eq $GENOTYPES_HEADER} split("\t", $header);
    $Genotypes_Col == -1 and die "$GENOTYPES_HEADER not found in header";
    
    while (<$fh>) {
        my @fields = split("\t");
        my $info = $fields[$INFO_COL];
        
        #convert Info field into HASH table
        my %info_hash;
        for (split(/;/, $info)) {
            my ($key, $value) = split("=", $_);
            $info_hash{$key} = $value;
        }

        # get the number of Het and Hom genotypes
        my @genotypes = split(",", $fields[$Genotypes_Col]);
        my @het_genos = grep {m"0/1" or m"1/0"} @genotypes;
        my @hom_genos = grep {m"1/1"} @genotypes; 
        my $allele_count = 2 * @hom_genos + @het_genos;
        s|\(\d/\d\)|| for (@genotypes);
        
        # get the Seattle Seq and Annovar Genes
        my (@SS_genes, @Annovar_genes); 
        if ($info_hash{'GL'} and $info_hash{'GL'} ne ".") { 
            @SS_genes = split(qr"[,/\\]", $info_hash{'GL'});
        }
        if ($info_hash{'Gene.refGene'} and $info_hash{'Gene.refGene'} ne ".") { 
            @Annovar_genes = split(qr"[,/\\]", $info_hash{'Gene.refGene'});
        } 
        push (my @genes, @SS_genes, @Annovar_genes);
        
        # get the CADD and Polyphen scores
        my ($CADD, $polyphen);
        if (looks_like_number($info_hash{'CADD'})) { $CADD = $info_hash{'CADD'} } 
        if (looks_like_number($info_hash{'PH'})) { $polyphen = $info_hash{'PH'} } 
                 
        # set Functional / Deleterious Flags
        my ($isTruncating, $isSplice, $isFrameshift, $isMissense, $isDeleterious) = (0, 0, 0, 0, 0);
        if ($info =~ /stop.?gain|stop.?los|start.?gain|start.?los/i) {
            $isTruncating = 1;
        }
        elsif ($info =~ /splice/i) {$isSplice = 1;}
        elsif ($info =~ /[^n]frame.?shift/i) {$isFrameshift = 1;}
        elsif ($info =~ /missense/i) {$isMissense = 1;}
        
#         SIFT_pred=T
#         LRT_pred=D
#         MutationTaster_pred=D
#         MutationAssessor_pred=L
#         FATHMM_pred=T
#         RadialSVM_pred=T
#         LR_pred=T

        # Polyphen D = deleterious; P = probably deleterious
        # zgrep -oP "Polyphen2_HDIV_pred=." 3001_2_26_28_30.SS_Ann_Eff.tsv.gz | sort | uniq
        if ($info_hash{'Polyphen2_HDIV_pred'} =~ /[DP]/ || $info_hash{'Polyphen2_HVAR_pred'} =~ /[DP]/) {
            $isDeleterious = 1;
        }
        
        # create a hash of gene variants data for this file and globally for all files
        for my $g (uniq @genes) {
            $file_genes{$g} = 1;
            $genes{$g}{'variants'}++;
            $genes{$g}{'alt_alleles'} += $allele_count;
            $genes{$g}{'het_genos'} += @het_genos;
            $genes{$g}{'hom_genos'} += @hom_genos;
            $genes{$g}{'deleterious'} += $isDeleterious;
            $genes{$g}{'truncating'} += $isTruncating;
            $genes{$g}{'splice'}     += $isSplice;
            $genes{$g}{'frameshift'} += $isFrameshift;
            $genes{$g}{'missense'}   += $isMissense;
            $genes{$g}{'families'}{$fam}{'truncating'} += $isTruncating;
            $genes{$g}{'families'}{$fam}{'splice'}     += $isSplice;
            $genes{$g}{'families'}{$fam}{'frameshift'} += $isFrameshift;
            $genes{$g}{'families'}{$fam}{'missense'}   += $isMissense;
            
            for my $id (@genotypes) {
                $genes{$g}{'individuals'}{$id} = 1;
                $Individuals{$id} = 1;
            }
            
            if (!$genes{$g}{'CADD'} or ($CADD && $CADD > $genes{$g}{'CADD'})) { 
                $genes{$g}{'CADD'} = $CADD
            }
            if (!$genes{$g}{'polyphen'} or ($polyphen && $polyphen > $genes{$g}{'polyphen'})) { 
                $genes{$g}{'polyphen'} = $polyphen
            }
        }
    }
    
    my $gene_count = keys %file_genes;
    say STDERR "$gene_count genes in $file";
    
    close $fh;
    $File_count++;
}

my $gene_count = keys %genes;
my $id_count = keys %Individuals;
say STDERR "$gene_count unique genes in $id_count individuals in $File_count files";

my @header = qw(
  Gene  Variants  Del_Var  Fams  Individuals  Alt_Alleles  Het  Hom_Alt  
  LOF  Trunc  Splice  Frame  Mis  Var_per_1kB  CADDmax  PhenMax  FamIDS  ==);
say join("\t", @header);

for my $g (sort keys %genes) {
    my $length = $gene_lengths{$g};
    my $var_per_bp = $length ? sprintf("%.2f", $genes{$g}{'alt_alleles'} * 1000 / $length) : $DATA_NA;
    my @individuals = keys $genes{$g}{'individuals'};
    my $lof = $genes{$g}{'truncating'} + $genes{$g}{'splice'} + $genes{$g}{'frameshift'};
    
    my @fams;
    for (keys $genes{$g}{'families'}) {
        my $hfam = $genes{$g}{'families'}{$_};
        my ($t, $s, $f, $m) = ($$hfam{'truncating'}, $$hfam{'splice'}, $$hfam{'frameshift'}, $$hfam{'missense'});
        my $family = "$_=T$t:S$s:F$f:M$m";
        push @fams, $family;
    }

    my @output = ($g, #gene
                  $genes{$g}{'variants'},
                  $genes{$g}{'deleterious'},
                  scalar @fams,         #number of fams
                  scalar @individuals,  #number of individuals
                  $genes{$g}{'alt_alleles'},
                  $genes{$g}{'het_genos'}, 
                  $genes{$g}{'hom_genos'},
                  $lof,
                  $genes{$g}{'truncating'},
                  $genes{$g}{'splice'},
                  $genes{$g}{'frameshift'},
                  $genes{$g}{'missense'},
                  $var_per_bp, 
                  $genes{$g}{'CADD'} ? $genes{$g}{'CADD'} : $DATA_NA,
                  $genes{$g}{'polyphen'} ? $genes{$g}{'polyphen'} : $DATA_NA,
                  join(";", @fams),
                  "=="
    );
    say join("\t", @output);
}
