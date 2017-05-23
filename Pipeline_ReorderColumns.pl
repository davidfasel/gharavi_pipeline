#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments


# desired column order
my @columns = qw(
  CHROM  POS  ID  avsnp147  REF  ALT  QUAL  FILTER  FORMAT  GENOTYPES  
  Sample_ID(GT)  Samples_With_Variant  Hets  Homs  AD_Pass  Missing 
  DetailsAnn
  MAF_max  MostDeleterious    
  CADD_ann  CADD_MSC  
  Polyphen2_HVAR_score PolyPhen_MSC  
  SIFT_score  SIFT_MSC
  CLINSIG  CLNDBN
  
  GENES  
    KidDisorder KidComment KidInheritence MouseGene MouseTerm OMIM_Disorder Emerge EmergeSNP  
    RVI  RVI%  HI  HI_imp  HI%  HI%_imp GDI  GDI_Phred  GDI_Damage   
    LOFRare.al  TruncRare.al  FrameRare.al  SpliceRare.al  MisRare.al.Poly>0.9
    LOF.al  Trunc.al  Frame.al  Splice.al 
    Exp_LOF.var  N_LOF.var Z_LOF  pLI
    
  ==Annovar  GeneAnn FuncAnn Exonic DetailsAnn  SIFT_score  SIFT_pred Polyphen2_HDIV_score 
    Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred 
    MutationTaster_score MutationTaster_pred MutationAssessor_score MutationAssessor_pred 
    FATHMM_score FATHMM_pred PROVEAN_score PROVEAN_pred VEST3_score CADD_raw CADD_ann 
    DANN_score fathmm-MKL_coding_score fathmm-MKL_coding_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred 
    integrated_fitCons_score integrated_confidence_value GERP++_RS phyloP7way_vertebrate phyloP20way_mammalian 
    phastCons7way_vertebrate phastCons20way_mammalian SiPhy_29way_logOdds
  
  ==snpEFF  GeneEFF  FuncEFF  DetailsEFF
  
  ==Freq All_ESP EUR_ESP AFR_ESP All_1KG Afr_1KG Amr_1KG Eas_1KG Eur_1KG Sas_1KG 
    ALL_Exac AFR_Exac AMR_Exac EAS_Exac FIN_Exac NFE_Exac OTH_Exac SAS_Exac  
    ALL_Exac_nontcga  AFR_Exac_nontcga  AMR_Exac_nontcga  EAS_Exac_nontcga  FIN_Exac_nontcga  NFE_Exac_nontcga  OTH_Exac_nontcga  SAS_Exac_nontcga
);
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

