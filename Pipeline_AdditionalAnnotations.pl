#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);
use Getopt::Long;  #process command line arguments

my $GENE_HEADER = "GENES";
my $DATA_NA = ".";
my $PUBLIC_DB = "/media/Data/public_databases";

my $KIDNEY_GENES = "$PUBLIC_DB/GeneLists/KidneyGenes_20170322.txt";
my $MOUSE_GENES  = "$PUBLIC_DB/GeneLists/CAKUT_genes-mouse.txt";
my $OMIM         = "$PUBLIC_DB/GeneLists/OMIM_genemap.txt";
my $HI           = "$PUBLIC_DB/HI/HI_prediction.bed";
my $HI_IMP       = "$PUBLIC_DB/HI/HI_prediction_with_imputation.bed";
my $RVI          = "$PUBLIC_DB/RVI/Genic_Intolerance_Scores.txt";
my $GDI          = "$PUBLIC_DB/GDI/GDI_full_10282015.txt";
my $EXAC         = "$PUBLIC_DB/Exac/exac_LOF_alleles_by_gene.txt";
my $CONSTRAINT   = "$PUBLIC_DB/Exac/forweb_cleaned_exac_r03_march16_z_data_pLI.txt";
my $EMERGE       = "$PUBLIC_DB/GeneLists/EmergeGenes.txt";
my $EMERGE_SNP   = "$PUBLIC_DB/GeneLists/EmergeSNPs.tsv";
my $CADD_MSC       = "$PUBLIC_DB/GeneLists/MSC_CADD_95.tsv";
my $SIFT_MSC       = "$PUBLIC_DB/GeneLists/MSC_Sift_95.tsv";
my $Poly_MSC       = "$PUBLIC_DB/GeneLists/MSC_PolyPhen_95.tsv";
my $PHARMO         = "$PUBLIC_DB/GeneLists/pharmocogenetics.tsv";


# IMPORTANT: update these header names in main Pipeline
my @header = qw( 
  RVI  RVI%  HI  HI_imp  HI%  HI%_imp
  GDI	GDI_Phred	GDI_Damage
  KidDisorder  KidComment  KidInheritence  
  MouseGene  MouseTerm
  OMIM_Disorder  
  LOFRare.al  TruncRare.al  FrameRare.al  SpliceRare.al  MisRare.al.Poly>0.9
  LOF.al  Trunc.al  Frame.al  Splice.al 
  Exp_LOF.var  N_LOF.var Z_LOF  pLI
  CADD_MSC  PolyPhen_MSC  SIFT_MSC
  Emerge  EmergeSNP
  Pharmo
);

#

my %Annotations;
for my $h (@header) {
    $Annotations{$h} = {};
}

my (%evsburden, %rvi);

# Hash Kidney disorders file by gene name
# 1636	ACE	angiotensin I converting enzyme (peptidyl-dipeptidase A) 1	17q23	Renal tubular dysgenesis	(+)106180 		recessive
open(FILE, $KIDNEY_GENES) or die "Unable to open Kidney File: $KIDNEY_GENES";
while(my $line = <FILE>){
    chomp($line);
    my @kidney  = split(/\t/, $line, -1);
    my $gene = $kidney[1];
    $gene =~ s/\s+//g;
    $Annotations{"KidDisorder"}{$gene} = $kidney[5];
    $Annotations{"KidInheritence"}{$gene} = $kidney[3];
    $Annotations{"KidComment"}{$gene}  = $kidney[6];
}
close FILE;

# Hash Mouse genes file
#Input  Input Type    	MGI Gene/Marker ID	Symbol	Name                              	Feature Type       	MP ID     	Term
#AHI1  	current symbol	MGI:87971         	Ahi1  	Abelson helper integration site 1 	protein coding gene	MP:0003675	kidney cysts
open(FILE, "$MOUSE_GENES") or die "Unable to open Mouse Genes: $MOUSE_GENES";
while(my $line = <FILE>){
    chomp($line);
    my @fields = split(/\t/, $line);
    my $gene = $fields[0];
    $gene =~ s/\s+//g;
    $Annotations{"MouseGene"}{$gene} = $fields[3];
    $Annotations{"MouseTerm"}{$gene} = $fields[7];
}
close FILE;

# Hash OMIM file by gene name.  Example of OMIM line:
#                      5-Gene             6 7                           8 9      10               11               12
# 1.270|10|7|10|1p34.3|GJB3, CX31, DFNA2B|C|Gap junction protein, beta-3||603324|REn, REc, Psh, A|same YAC as GJA4| |
#   13                                                                   14                                                       15
#   Erythrokeratodermia variabilis et progressiva, 133200 (3); Deafness,|Deafness, autosomal dominant, with peripheral neuropathy|GJB2/GJB3||
open(FILE,"$OMIM") or die "Unable to open OMIM file: $OMIM";
while(my $line = <FILE>){
    chomp($line);
    my ($genes, $comment, $disorder1, $disorder2, $disorder3) = ( split(/\|/, $line) )[5, 11, 13, 14, 15];
    my @gene_list = split(/, +/, $genes);   #gene names are separated by comma and space
    for my $gene (@gene_list){
        $disorder1 =~ s/^\s+$//g;
        my @out = grep {/\S+/} ($disorder1, $disorder2, $disorder3, $comment);  # get fields that aren't just whitespace
        $Annotations{"OMIM_Disorder"}{$gene} = join("|", @out) if $disorder1;
    }
}
close FILE;

open(FILE, "$RVI") or die "Unable to open RVI file: $RVI";
while(my $line = <FILE>){
    next if $. == 1;
    $line =~ s/\r?\n//;
    my ($gene, $rvi_score, $rvi_percent) = split(/\t/, $line);
    $Annotations{"RVI"}{$gene} = $rvi_score;
    $Annotations{"RVI%"}{$gene} = $rvi_percent;
}
close FILE;

# Hash HI files by gene name
#chr1	850392	869824	SAMD11|0.085|79.5%	0.085
open(FILE, "$HI") or die "Unable to open HI file: $HI";
while(my $line = <FILE>){
    next if $. == 1;  #skip track name line
    chomp $line;
    my ($gene, $hi_score, $hi_percent) = split(/\|/, (split /\t/, $line)[3]);
    $Annotations{"HI"}{$gene} = $hi_score;
    $Annotations{"HI%"}{$gene} = $hi_percent;
}
close FILE;
open(FILE, "$HI_IMP") or die "Unable to open HI_imp file: $HI_IMP";
while(my $line = <FILE>){
    next if $. == 1;  #skip track name line
    $line =~ s/\r?\n//;
    my ($gene, $hii_score, $hii_percent) = split(/\|/, (split /\t/, $line)[3]);
    $Annotations{"HI_imp"}{$gene} = $hii_score;
    $Annotations{"HI%_imp"}{$gene} = $hii_percent;
}
close FILE;

# Hash GDI file by gene name
#Gene	GDI	GDI-Phred	Gene damage prediction (all disease-causing genes
#DEFB104B	0.00011	0.00000	Low
open(FILE, "$GDI") or die "Unable to open GDI file: $GDI";
while(my $line = <FILE>){
    next if $. == 1;  #skip header
    $line =~ s/\r?\n//;
    my ($gene, $gdi_score, $gdi_phred, $gdi_pred) = split(/\t/, $line);
    $Annotations{"GDI"}{$gene}       = $gdi_score;
    $Annotations{"GDI_Phred"}{$gene} = $gdi_phred;
    $Annotations{"GDI_Damage"}{$gene} = $gdi_pred;
}
close FILE;

# My summary file of Exac alleles by gene created by exac_allele_counts.pl
# Gene    LOF     Trunc   Frame   Splice  MisPoly>0.9
open(FILE, "$EXAC") or die "Unable to open Exac Alleles file: $EXAC";
while(my $line = <FILE>){
    next if $. == 1;  #skip header
    $line =~ s/\r?\n//;
    my ($gene, $lof_r, $trunc_r, $frame_r, $splice_r, $mis, $lof, $trunc, $frame, $splice) = split(/\t/, $line);
    $Annotations{"LOFRare.al"}{$gene}    = $lof_r;
    $Annotations{"TruncRare.al"}{$gene}  = $trunc_r;
    $Annotations{"FrameRare.al"}{$gene}  = $frame_r;
    $Annotations{"SpliceRare.al"}{$gene} = $splice_r;
    $Annotations{"MisRare.al.Poly>0.9"}{$gene} = $mis;
    $Annotations{"LOF.al"}{$gene}        = $lof;
    $Annotations{"Trunc.al"}{$gene}      = $trunc;
    $Annotations{"Frame.al"}{$gene}      = $frame;
    $Annotations{"Splice.al"}{$gene}     = $splice;
}
close FILE;

# My summary file of Exac alleles by gene created by exac_allele_counts.pl
#0         	1   	2  	3      	4       	5     	6 	7     	8     	9     	10   	11   	12   	13     	14     	15     	16   	17   	18   	19
#transcript	gene	chr	n_exons	tx_start	tx_end	bp	mu_syn	mu_mis	mu_lof	n_syn	n_mis	n_lof	exp_syn	exp_mis	exp_lof	syn_z	mis_z	lof_z	pLI
open(FILE, "$CONSTRAINT") or die "Unable to open Exac Constraint file: $CONSTRAINT";
while(my $line = <FILE>){
    next if $. == 1;  #skip header
    $line =~ s/\r?\n//;
    my @fields = split(/\t/, $line);
    my ($gene, $exp_lof, $n_lof, $lof_z, $pli) = @fields[1, 15, 12, 18, 19];
    $Annotations{"Exp_LOF.var"}{$gene}   = $exp_lof;
    $Annotations{"N_LOF.var"}{$gene}     = $n_lof;
    $Annotations{"Z_LOF"}{$gene}         = $lof_z;
    $Annotations{"pLI"}{$gene}           = $pli;
}
close FILE;

open(FILE, "$EMERGE") or die "Unable to open Emerge file: $EMERGE";
while(my $line = <FILE>){
    $line =~ s/\r?\n//;
    my ($gene, $association) = split(/\t/, $line);
    $Annotations{"Emerge"}{$gene} = $association;
}
close FILE;

open(FILE, "$EMERGE_SNP") or die "Unable to open Emerge file: $EMERGE_SNP";
while(my $line = <FILE>){
    chomp $line;
    my ($rsID, $location, $source) = split(/\t/, $line);
    $Annotations{"EmergeSNP"}{$location} = "$rsID,$source";
}
close FILE;

open(FILE, "$CADD_MSC") or die "Unable to open CADD_MSC file: $CADD_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $Annotations{"CADD_MSC"}{$gene} = $msc;
}
close FILE;

open(FILE, "$Poly_MSC") or die "Unable to open PolyPhen_MSC file: $Poly_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $Annotations{"PolyPhen_MSC"}{$gene} = $msc;
}
close FILE;

open(FILE, "$SIFT_MSC") or die "Unable to open SIFT_MSC file: $SIFT_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $Annotations{"SIFT_MSC"}{$gene} = $msc;
}
close FILE;

open(FILE, "$PHARMO") or die "Unable to open PHARMO file: $PHARMO";
while(my $line = <FILE>){
    chomp $line;
    
    my ($drug, $area, $gene_clean, $gene, $type) = split(/\t/, $line);
    
    my @gene_list = split(/,/, $gene_clean);   #gene names are separated by comma and possible space
    s/ +// for @gene_list;   #remove white space around gene names
    for my $gene (@gene_list){
        my $out = join("|", ($drug, $area, $gene, $type));
        # since there could be more than one drug for a gene we need to append the fields
        $Annotations{"Pharmo"}{$gene} .= "$out;  ";
    }
}
close FILE;

for my $key (keys %Annotations) {
    $Annotations{$key}{"output"} = [];
}



#print the header
my $line = <STDIN>;
chomp $line;
my @fields = split(/\t/, $line);
my @output = ($line, @header);
say join("\t", @output);

# get index of Gene Column
my $Gene_Col = first_index {$_ eq $GENE_HEADER} @fields ;
die "unable to find Genes Column" if ($Gene_Col == -1);


### print each line with the new annotations
while(my $line = <STDIN>) {
    chomp($line);
    my @output = ($line);
    my @fields = split(/\t/,$line);
      
    #get genes from GENES column
    my @genes = split(",", $fields[$Gene_Col]);
    s/"//g for @genes;  # the gene names my be surrounded by quotes as an artifact from Excel

    # get the location of the variant from the first 2 columns
    my ($chrom, $position) = @fields[0,1];
    $chrom =~ s/chr//i;  #remove chr prefix if present
    my $loc = "chr$chrom:$position";
    
    for my $gene (@genes){
        for my $key (keys %Annotations) {
            if(exists $Annotations{$key}{$gene}){
                push (@{$Annotations{$key}{"output"}}, $Annotations{$key}{$gene});
            }
        }
    }

    # going through each header we get the annotation for that column
    for my $h (@header) {
        @{$Annotations{$h}{"output"}} or @{$Annotations{$h}{"output"}} = ($DATA_NA);
        push( @output, join(";", @{$Annotations{$h}{"output"}}) );
    }
    
    # check for Emerge pathogenic SNPs
    $Annotations{"EmergeSNP"} or $Annotations{"EmergeSNP"} = ".";
    push( @output, $Annotations{"EmergeSNP"});
    
    say join("\t", @output);
    
    # clear the output arrays to be read for the next line
    for my $key (keys %Annotations) {
        $Annotations{$key}{"output"} = [];
    }
    
}
close FILE;
