#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use List::MoreUtils qw(first_index);
use Getopt::Long;  #process command line arguments

my $GENE_HEADER = "GENES";
my $DATA_NA = ".";
my $PUBLIC_DB = "/media/Data/public_databases";

my $KIDNEY_GENES = "$PUBLIC_DB/GeneLists/KidneyGenes6.txt";
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



# IMPORTANT: update these header names in main Pipeline
my @header = qw( 
  RVI  RVI%  HI  HI_imp  HI%  HI%_imp
  GDI	GDI_Phred	GDI_Damage
  KidDisorder  KidComment  KidInheritence  
  MouseGene  MouseTerm
  OMIM_Disorder  Emerge
  LOFRare.al  TruncRare.al  FrameRare.al  SpliceRare.al  MisRare.al.Poly>0.9
  LOF.al  Trunc.al  Frame.al  Splice.al 
  Exp_LOF.var  N_LOF.var Z_LOF  pLI
  CADD_MSC  PolyPhen_MSC  SIFT_MSC

);

my %hash_all;
for my $h (@header) {
    $hash_all{$h} = {};
}



my (%kidney_disorder, %mouse_genes, %omim, %evsburden, %rvi, %hi, %hii, 
      %gdi, %exac, %constraint, %emerge, %msc_sift, %msc_poly);

# Hash Kidney disorders file by gene name
# 1636	ACE	angiotensin I converting enzyme (peptidyl-dipeptidase A) 1	17q23	Renal tubular dysgenesis	(+)106180 		recessive
open(FILE, $KIDNEY_GENES) or die "Unable to open Kidney File: $KIDNEY_GENES";
while(my $line = <FILE>){
    chomp($line);
    my $gene = (split(/\t/,$line))[1];
    $gene =~ s/\s+//g;
    $kidney_disorder{$gene} = $line;
}
close FILE;

# Hash Mouse genes file
#Input  Input Type    	MGI Gene/Marker ID	Symbol	Name                              	Feature Type       	MP ID     	Term
#AHI1  	current symbol	MGI:87971         	Ahi1  	Abelson helper integration site 1 	protein coding gene	MP:0003675	kidney cysts
open(FILE, "$MOUSE_GENES") or die "Unable to open Mouse Genes: $MOUSE_GENES";
while(my $line = <FILE>){
    chomp($line);
    my @fields = split(/\t/, $line, -1);
    my $gene = $fields[0];
    $gene =~ s/\s+//g;
    $mouse_genes{$gene} = [@fields[3,7]];
}
close FILE;

# Hash OMIM file by gene name
# 1.1|9|11|95|1pter-p36.13|CCV|P|Cataract, congenital, Volkmann type||115665|Fd|linked to Rh in Scottish family||Cataract, congenital, Volkmann type (2)| | ||
open(FILE,"$OMIM") or die "Unable to open OMIM file: $OMIM";
while(my $line = <FILE>){
    chomp($line);
    my @genes_omim = split(/\,/, (split(/\|/,$line))[5]);
    for my $gene (@genes_omim){
        $gene =~ s/\s+//g;
        $omim{$gene} = $line;
    }
}
close FILE;

# Hash RVI file by gene name
open(FILE, "$RVI") or die "Unable to open RVI file: $RVI";
while(my $line = <FILE>){
    next if $. == 1;
    $line =~ s/\r?\n//;
    my ($gene, $rvi_score, $rvi_percent) = split(/\t/, $line);
    $rvi{$gene} = [$rvi_score, $rvi_percent];
}
close FILE;


# Hash HI file by gene name
#chr1	850392	869824	SAMD11|0.085|79.5%	0.085
open(FILE, "$HI") or die "Unable to open HI file: $HI";
while(my $line = <FILE>){
    next if $. == 1;  #skip track name line
    $line =~ s/\r?\n//;
    my ($gene, $hi_score, $hi_percent) = split(/\|/, (split /\t/, $line)[3]);
    $hi{$gene} = [$hi_score, $hi_percent];
}
close FILE;
open(FILE, "$HI_IMP") or die "Unable to open HI_imp file: $HI_IMP";
while(my $line = <FILE>){
    next if $. == 1;  #skip track name line
    $line =~ s/\r?\n//;
    my ($gene, $hii_score, $hii_percent) = split(/\|/, (split /\t/, $line)[3]);
    $hii{$gene} = [$hii_score, $hii_percent];
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
    $gdi{$gene} = [$gdi_score, $gdi_phred, $gdi_pred];
}
close FILE;

# My summary file of Exac alleles by gene created by exac_allele_counts.pl
#Gene    LOF     Trunc   Frame   Splice  MisPoly>0.9
open(FILE, "$EXAC") or die "Unable to open Exac Alleles file: $EXAC";
while(my $line = <FILE>){
    next if $. == 1;  #skip header
    $line =~ s/\r?\n//;
    my ($gene, $lof_r, $trunc_r, $frame_r, $splice_r, $mis, $lof, $trunc, $frame, $splice) = split(/\t/, $line);
    $exac{$gene} = [$lof_r, $trunc_r, $frame_r, $splice_r, $mis, $lof, $trunc, $frame, $splice];
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
    $constraint{$gene} = [$exp_lof, $n_lof, $lof_z, $pli];
}
close FILE;

open(FILE, "$EMERGE") or die "Unable to open EMERGE file: $EMERGE";
while(my $line = <FILE>){
    $line =~ s/\r?\n//;
    my ($gene, $association) = split(/\t/, $line);
    $emerge{$gene} = $association;
}
close FILE;

open(FILE, "$CADD_MSC") or die "Unable to open CADD_MSC file: $CADD_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $hash_all{"CADD_MSC"}{$gene} = $msc;
    $hash_all{"CADD_MSC"}{"output"} = [];
}
close FILE;

# open(FILE, "$Poly_MSC") or die "Unable to open PolyPhen_MSC file: $Poly_MSC";
# while(my $line = <FILE>){
#     chomp $line;
#     my ($gene, $msc) = (split(/\t/, $line))[0,8];
#     $hash_all{"PolyPhen_MSC"}{$gene} = $msc;
#     $hash_all{"PolyPhen_MSC"}{"output"} = [];
# }
# close FILE;

open(FILE, "$Poly_MSC") or die "Unable to open PolyPhen_MSC file: $Poly_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $msc_poly{$gene} = $msc;
}
close FILE;

open(FILE, "$SIFT_MSC") or die "Unable to open SIFT_MSC file: $SIFT_MSC";
while(my $line = <FILE>){
    chomp $line;
    my ($gene, $msc) = (split(/\t/, $line))[0,8];
    $msc_sift{$gene} = $msc;
}
close FILE;

# for my $key (keys %hash_all) {
#     $hash_all{$key}{"output"} = [];
# }



#print the header
my $line = <STDIN>;
chomp $line;
my @fields = split(/\t/, $line);
my @output = ($line, @header);
say join("\t", @output);

# get index of Gene Column
my $Gene_Col = first_index {$_ eq $GENE_HEADER} @fields ;
die "unable to find Genes Column" if ($Gene_Col == -1);
  
while(my $line = <STDIN>) {
    chomp($line);
    my @output = ($line);
    my @fields = split(/\t/,$line);
      
    #get genes from GENES column
    my @genes = split(",", $fields[$Gene_Col]);
    s/"//g for @genes;  # the gene names my be surrounded by quotes as an artifact from Excel

    #goto next line if no genes
    if (not @genes) {
        for my $h (@header) {
            push (@output, $h eq "==" ? "==" : $DATA_NA );
        }
        say join("\t", @output);
        next;
    }
    
    # get the annotations for each gene (since there could be more than one gene)
    my (@rvi_scores, @rvi_percentiles, @hi_scores, @hii_scores, @hi_percentiles, @hii_percentiles, 
        @gdi_scores, @gdi_phreds, @gdi_preds,
        @kid_disorder, @kid_comment, @kid_inher, 
        @mouse_gene, @mouse_term,  
        @omim_disorders,  @emerge_ass,
        @exac_lof_r, @exac_trunc_r, @exac_frame_r, @exac_splice_r, @exac_mis,
        @exac_lof, @exac_trunc, @exac_frame, @exac_splice, 
        @const_explof, @const_nlof, @const_lofz, @const_pli, 
        @msc_cadds, @msc_polys, @msc_sifts);
      
#     my @all_arrays = (
#         \@rvi_scores, \@rvi_percentiles, \@hi_scores, \@hii_scores, \@hi_percentiles, \@hii_percentiles, 
#         \@gdi_scores, \@gdi_phreds, \@gdi_preds,
#         \@kid_disorder, \@kid_comment, \@kid_inher, 
#         \@mouse_gene, \@mouse_term,  
#         \@omim_disorders, \@emerge_ass,
#         \@exac_lof_r, \@exac_trunc_r, \@exac_frame_r, \@exac_splice_r, \@exac_mis,
#         \@exac_lof, \@exac_trunc, \@exac_frame, \@exac_splice, 
#         \@const_explof, \@const_nlof, \@const_lofz, \@const_pli, 
#         \@msc_cadds, \@msc_polys, \@msc_sifts
#     );
    
    my @all_arrays = (
        \@rvi_scores, \@rvi_percentiles, \@hi_scores, \@hii_scores, \@hi_percentiles, \@hii_percentiles, 
        \@gdi_scores, \@gdi_phreds, \@gdi_preds,
        \@kid_disorder, \@kid_comment, \@kid_inher, 
        \@mouse_gene, \@mouse_term,  
        \@omim_disorders, \@emerge_ass,
        \@exac_lof_r, \@exac_trunc_r, \@exac_frame_r, \@exac_splice_r, \@exac_mis,
        \@exac_lof, \@exac_trunc, \@exac_frame, \@exac_splice, 
        \@const_explof, \@const_nlof, \@const_lofz, \@const_pli, 
        $hash_all{"CADD_MSC"}{"output"}, \@msc_polys, \@msc_sifts
    );
    
#    $hash_all{"Poly_MSC"}{"output"}
    
    for my $gene (@genes){
        if(exists $kidney_disorder{$gene}){
            my @kidney = split(/\t/, $kidney_disorder{$gene});
            push (@kid_disorder, $kidney[4]);
            push (@kid_inher, $kidney[7]);
            push (@kid_comment, $kidney[8]);
        }
        if(exists $mouse_genes{$gene}) {
            push(@mouse_gene, $mouse_genes{$gene}[0]);
            push(@mouse_term, $mouse_genes{$gene}[1]);
        }
        if(exists $omim{$gene}) {
            my @omim_val = split(/\|/, $omim{$gene});
            push(@omim_disorders, $omim_val[13], $omim_val[14], $omim_val[15]);
        }
        if(exists $rvi{$gene}) {
            push (@rvi_scores, $rvi{$gene}[0]);
            push (@rvi_percentiles, $rvi{$gene}[1]);
        }
        if(exists $hi{$gene}) {
            push (@hi_scores, $hi{$gene}[0]);
            push (@hi_percentiles, $hi{$gene}[1]);
        }
        if(exists $hii{$gene}) {
            push (@hii_scores, $hii{$gene}[0]);
            push (@hii_percentiles, $hii{$gene}[1]);
        }
        if(exists $gdi{$gene}) {
            push (@gdi_scores, $gdi{$gene}[0]);
            push (@gdi_phreds, $gdi{$gene}[1]);
            push (@gdi_preds, $gdi{$gene}[2]);
        }
        if(exists $exac{$gene}) {
            push (@exac_lof_r, $exac{$gene}[0]);
            push (@exac_trunc_r, $exac{$gene}[1]);
            push (@exac_frame_r, $exac{$gene}[2]);
            push (@exac_splice_r, $exac{$gene}[3]);
            push (@exac_mis, $exac{$gene}[4]);
            push (@exac_lof, $exac{$gene}[5]);
            push (@exac_trunc, $exac{$gene}[6]);
            push (@exac_frame, $exac{$gene}[7]);
            push (@exac_splice, $exac{$gene}[8]);
        }
        if(exists $constraint{$gene}) {
            push (@const_explof, $constraint{$gene}[0]);
            push (@const_nlof, $constraint{$gene}[1]);
            push (@const_lofz, $constraint{$gene}[2]);
            push (@const_pli, $constraint{$gene}[3]);
        } 
        if(exists $emerge{$gene}) {
            push (@emerge_ass, $emerge{$gene});
        } 
        if(exists $hash_all{"CADD_MSC"}{$gene}) {
            push (@{$hash_all{"CADD_MSC"}{"output"}}, $hash_all{"CADD_MSC"}{$gene});
        }
#         if(exists $hash_all{"Poly_MSC"}{$gene}) {
#             push (@{$hash_all{"Poly_MSC"}{"output"}}, $hash_all{"Poly_MSC"}{$gene});
#         }
       if(exists $msc_poly{$gene}) {
            push (@msc_polys, $msc_poly{$gene});
        }
        if(exists $msc_sift{$gene}) {
            push (@msc_sifts, $msc_sift{$gene});
        }
    }
    #remove any that are just whitespace
    @kid_disorder   = grep { /\S/ } @kid_disorder;
    @kid_inher      = grep { /\S/ } @kid_inher;
    @kid_comment    = grep { /\S/ } @kid_comment;
    @mouse_gene     = grep { /\S/ } @mouse_gene;
    @mouse_term     = grep { /\S/ } @mouse_term;
    @omim_disorders = grep { /\S/ } @omim_disorders;
    
    for (@all_arrays) {
        @{$_} or @{$_} = ($DATA_NA);   #if the array is empty, set it to the Missing value
        push(@output, join(";", @{$_}));
    }
    
    say join("\t", @output);
    
    # clear the output arrays to be read for the next line
    for my $key (keys %hash_all) {
        $hash_all{$key}{"output"} = [];
    }
    
    
    
}
close FILE;
