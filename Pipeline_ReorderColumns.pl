#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments


# desired column order
my @columns = qw(
  CHROM  POS avsnp150  REF  ALT  QUAL  FILTER  
  Sample_ID(GT)  Samples_With_Variant  
  MAF_max  MostDeleterious 
  
  GENES  
    KidDisorder KidComment KidInheritence MouseGene MouseTerm OMIM_Disorder 
    Emerge EmergeSNP  Pharma
    
  FORMAT  GENOTYPES Hets  Homs  AD_Pass  Missing 

  ==snpEFF 
    snpEffGene Annotation Impact LOF

  ==Annovar  
    AAChange.refGene 
    Func.refGene 
    Gene.refGene 
    GeneDetail.refGene 
    ExonicFunc.refGene 
    CADD_phred 
    CADD_raw 
    CADD_raw_rankscore 
    CLNALLELEID 
    CLNDISDB 
    CLNDN 
    CLNREVSTAT 
    CLNSIG 
    DANN_rankscore 
    DANN_score 
    Eigen_coding_or_noncoding 
    Eigen-PC-raw 
    Eigen-raw 
    FATHMM_converted_rankscore 
    FATHMM_pred 
    FATHMM_score 
    fathmm-MKL_coding_pred 
    fathmm-MKL_coding_rankscore 
    fathmm-MKL_coding_score 
    GenoCanyon_score 
    GenoCanyon_score_rankscore 
    GERP++_RS 
    GERP++_RS_rankscore 
    GTEx_V6p_gene 
    GTEx_V6p_tissue 
    integrated_confidence_value 
    integrated_fitCons_score 
    integrated_fitCons_score_rankscore 
    Interpro_domain 
    LRT_converted_rankscore 
    LRT_pred 
    LRT_score 
    M-CAP_pred 
    M-CAP_rankscore 
    M-CAP_score 
    MetaLR_pred 
    MetaLR_rankscore 
    MetaLR_score 
    MetaSVM_pred 
    MetaSVM_rankscore 
    MetaSVM_score 
    MutationAssessor_pred 
    MutationAssessor_score 
    MutationAssessor_score_rankscore 
    MutationTaster_converted_rankscore 
    MutationTaster_pred 
    MutationTaster_score 
    MutPred_rankscore 
    MutPred_score 
    phastCons100way_vertebrate 
    phastCons100way_vertebrate_rankscore 
    phastCons20way_mammalian 
    phastCons20way_mammalian_rankscore 
    phyloP100way_vertebrate 
    phyloP100way_vertebrate_rankscore 
    phyloP20way_mammalian 
    phyloP20way_mammalian_rankscore 
    Polyphen2_HDIV_pred 
    Polyphen2_HDIV_rankscore 
    Polyphen2_HDIV_score 
    Polyphen2_HVAR_pred 
    Polyphen2_HVAR_rankscore 
    Polyphen2_HVAR_score 
    PROVEAN_converted_rankscore 
    PROVEAN_pred 
    PROVEAN_score 
    REVEL_rankscore 
    REVEL_score 
    SIFT_converted_rankscore 
    SIFT_pred 
    SIFT_score 
    SiPhy_29way_logOdds 
    SiPhy_29way_logOdds_rankscore 
    VEST3_rankscore 
    VEST3_score 
  
  ==MAFs  
    G_AF G_AF_raw G_AF_male G_AF_female G_AF_afr G_AF_ami G_AF_amr G_AF_asj G_AF_eas G_AF_fin G_AF_nfe G_AF_oth G_AF_sas 
    GME_NWA  GME_NEA  GME_AP  GME_Israel  GME_SD  GME_TP  GME_CA

  ==GeneScores
  GENES
  RVI  RVI%  HI  HI_imp  HI%  HI%_imp GDI  GDI_Phred  GDI_Damage   
    LOFRare.al  TruncRare.al  FrameRare.al  SpliceRare.al  MisRare.al.Poly>0.9
    LOF.al  Trunc.al  Frame.al  Splice.al 
    Exp_LOF.var  N_LOF.var Z_LOF  pLI
);

# the following are included in Gnomad
# All_ESP EUR_ESP AFR_ESP All_1KG Afr_1KG Amr_1KG Eas_1KG Eur_1KG Sas_1KG 

#SVM_PROBABILITY  SVM_POSTERIOR

#   ==SeattleSeq  geneList  functionGVS  functionDBSNP  accession  aminoAcids  proteinPosition
#     cDNAPosition  chimpAllele  clinicalAssociation  distanceToSplice  keggPathway  tfbs
#     PPI  proteinAccession  CADD  PolyPhen  PhastCons  GERPConsScore  granthamScore  microRNAs


my $line = <STDIN>;
chomp $line;
my @headers = split(/\t/, $line);
my @indexes;


# get the index of each of the columns based on the header names above
# put these indexes into an array in the order that they appear above
for my $col (@columns) {
  my $index = 0;
  ++$index until $headers[$index] eq $col or $index > $#headers;
  #first_index{$col eq $_} @headers;
  die "Can't reorder columns, couldn't find column header $col" if $index > $#headers;
  push(@indexes, $index);
}


# print the header
say join("\t", @columns);

while (<STDIN>) {
  chomp;
  my @fields = split("\t");
  my @out;
  push(@out, $fields[$_]) for (@indexes);  
  say join("\t", @out);  
}

