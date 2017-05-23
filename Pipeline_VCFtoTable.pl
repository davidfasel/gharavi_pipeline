#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my $MISSING = ".";  #define what should be shown if field is missing (can be empty string if desired)
my $ALLELIC_DEPTH_RATIO = 0.3;

my $MAF_filter = 1;
my $vcf_file;
GetOptions (
    'file|f=s' => \$vcf_file,
    'maf|m:f' => \$MAF_filter,
) or die "Usage: Pipeline_VCFtoTable.pl -m 0.01 -f file.vcf.gz";

my (%dbsnp);

# NOTE!  if any of these headers change, they need to be updated in the main pipeline file
my @header = qw(
    CHROM  POS  ID  avsnp147  REF  ALT  QUAL  FILTER  INFO  FORMAT GENOTYPES 
    Sample_ID(GT) Samples_With_Variant  Hets  Homs  AD_Pass  Missing  MAF_max  
    
    SVM_PROBABILITY SVM_POSTERIOR
    
    GENES  CLINSIG  CLNDBN

    ==SeattleSeq
    geneList
    functionGVS
    functionDBSNP
    accession
    aminoAcids
    proteinPosition
    cDNAPosition
    chimpAllele
    clinicalAssociation
    distanceToSplice
    keggPathway
    tfbs
    PPI
    proteinAccession
    CADD  
    PolyPhen
    PhastCons  
    GERPConsScore
    granthamScore
    microRNAs
    
    ==Annovar  GeneAnn  FuncAnn  Exonic  DetailsAnn  
    SIFT_score  SIFT_pred
    Polyphen2_HDIV_score  Polyphen2_HDIV_pred  Polyphen2_HVAR_score  Polyphen2_HVAR_pred
    LRT_score  LRT_pred
    MutationTaster_score  MutationTaster_pred  MutationAssessor_score  MutationAssessor_pred
    FATHMM_score  FATHMM_pred
    PROVEAN_score  PROVEAN_pred
    VEST3_score
    CADD_raw  CADD_ann
    DANN_score
    fathmm-MKL_coding_score  fathmm-MKL_coding_pred
    MetaSVM_score  MetaSVM_pred  MetaLR_score  MetaLR_pred
    integrated_fitCons_score  integrated_confidence_value
    GERP++_RS
    phyloP7way_vertebrate  phyloP20way_mammalian
    phastCons7way_vertebrate  phastCons20way_mammalian
    SiPhy_29way_logOdds
    
    ==snpEFF  GeneEFF  FuncEFF  DetailsEFF
    
    ==Freq  All_ESP  EUR_ESP  AFR_ESP
    All_1KG  Afr_1KG  Amr_1KG  Eas_1KG  Eur_1KG  Sas_1KG
    ALL_Exac  AFR_Exac  AMR_Exac  EAS_Exac  FIN_Exac  NFE_Exac  OTH_Exac  SAS_Exac 
    ALL_Exac_nontcga  AFR_Exac_nontcga  AMR_Exac_nontcga  EAS_Exac_nontcga  FIN_Exac_nontcga  NFE_Exac_nontcga  OTH_Exac_nontcga  SAS_Exac_nontcga
    
   
);
                  
#    ==dbSNP  141_ID  141_REF  141_ALT  141_AF  AlleleStatus  
say join("\t", @header);
 


#### Build the Output file
if ($vcf_file =~ /\.gz$/) {
    open(FILE, '-|', 'gzip', '-dc', $vcf_file) or die "File not found: $vcf_file";
}
else {
    open(FILE, $vcf_file) or die "File not found: $vcf_file";
}



my @sample_ids;
while(my $line = <FILE>)
{   
    print STDERR "          Processing Line $. \r" if ($. % 1000 == 0);  
    
    chomp($line);
    my @var=split(/\t/,$line);
    
    # get the sample IDs
    if($line =~ m/^#CHROM/){
        for(my $s = 9; $s < @var; $s++){
            push (@sample_ids, $var[$s]);
        }
    }
    #skip comment or blank lines
    next if ($line =~ /^\s*$/ or $line =~ /^#/ );
    
    #convert Info field into HASH table
    my %info_hash;
    my @info_field = split(/;/, $var[7]);
    for my $item (@info_field) {
      my ($key, $value) = split("=", $item);
      #special case where item has no value, in which case it's a flag
      $value //= "flag";
      $info_hash{$key} = $value;
    }
    
    
    #### to filter by a different population, replace "PopFreqMax" below
    my $MAF_max = $info_hash{"PopFreqMax"} || $MISSING;
    next if looks_like_number($MAF_max) && $MAF_max > $MAF_filter; 

    # get the index of the depth (DP) value in the FORMAT field.  
    # Format field may not exist if this is just a reference of variants without genotypes (such as EXAC)
    my ($dp_index, $ad_index);
    my $format = ($var[8] or "");
    if ($format) {
        my @dp_col = split(/\:/,$var[8]);
        for(0 .. @dp_col-1) {
            $dp_index = $_ if ($dp_col[$_] eq "DP");
            $ad_index = $_ if ($dp_col[$_] eq "AD");
        }
    }

    # count the samples and get their ID's, Genotypes, and Depths, Number that are missing and 
    # whether Allelic Depth ratio is greater than ALLELIC_DEPTH_RATIO for hets
    my ($samp) = ("");
    my ($num_samples, $hets, $homs) = (0, 0, 0);
    my $genotypes = "";
    my $missing = 0;
    my $PASS_AD = 0;
    
    for (my $c = 9; $c < @var; $c++) {
        my @gt_fields = split(/\:/, $var[$c]);
        
        #split genotype by / or | (phased genotype)
        my @gt = split("[/\|]", $gt_fields[0]);  
                
        if ($gt[0] eq ".") {  
            $missing++;  # count missing genotypes
        }
        elsif ($gt_fields[0] !~ "^0.0" ) {
            $num_samples++;
            $gt[0] eq $gt[1] ? $homs++ : $hets++;
            
            $genotypes .= "$var[$c];";
            $samp .= $sample_ids[$c-9] . "($gt_fields[0]),"; 
            
            #check Allelic Depth ratio if heterozygous
            if ($ad_index && $gt_fields[$ad_index] ne ".") { 
                if ($gt[0] ne $gt[1]) {  #check for het
                    my @AD = split(",", $gt_fields[$ad_index]);
                    my $sum = $AD[0] + $AD[1];
                    $PASS_AD++ if ($sum > 0 && ($AD[1] / $sum) >= $ALLELIC_DEPTH_RATIO);
                }
                else {  # count as pass if homozygous
                  $PASS_AD++;
                }
            }

            # get the depths, prefer Allelic depth (AD) over just regular depth (DP)
#             if ($ad_index && $gt_fields[$ad_index] ne ".") {
#                 $depth .= "$gt_fields[$ad_index];";               
#             }
#             elsif ($dp_index && $gt_fields[$dp_index] ne ".") {
#                 $depth .= "$gt_fields[$dp_index];";
#             }
#             else {
#                 $depth .= ".;";
#             }

        }
    }
    next if $num_samples == 0;
    $genotypes = $samp = "> 100" if ($num_samples > 100);
    
    # get a list of all genes  todo:SnpEff
    my @genes;
    my $seaseq_genes = $info_hash{"GL"} || $MISSING; 
    my $annovar_genes = $info_hash{"Gene.refGene"} || $MISSING;    
        
    push(@genes, split("[,/]", $seaseq_genes), 
                 split(/,|\\x3b/, $annovar_genes));
    @genes = uniq(grep {$_ ne "."} @genes);
    my $gene = join(",", @genes) . ",";  # append a comma to prevent excel from converting some genes to dates
    
    my $clinsig = $info_hash{"CLINSIG"} || $MISSING;
    my $clindbn = $info_hash{"CLNDBN"} || $MISSING;
    my $clnacc = $info_hash{"CLNACC"} || $MISSING;
    my $clndsdb = $info_hash{"CLNDSDB"} || $MISSING;
    my $clndsdbid = $info_hash{"CLNDSDBID"} || $MISSING;
    my $clinout = join("|", $clindbn, $clnacc, $clndsdb, $clndsdbid);
    $clinout =~ s/\\x2c/,/g;  #commas are encoded, so change back to comma
    $clinout = "." if $clinsig eq ".";
    
    my @out = (
      $var[0], $var[1],  $var[2],
      $info_hash{"avsnp147"} || $MISSING,
      $var[3], $var[4], $var[5], $var[6], $var[7],
      $format,
      $genotypes,
      $samp,
      $num_samples,
      $hets,
      $homs,
      $PASS_AD,
      $missing,
      $MAF_max,

      $info_hash{'SVM_PROBABILITY'} || $MISSING, 
      $info_hash{'SVM_POSTERIOR'} || $MISSING,
      
      $gene,
      $clinsig,
      $clinout,
    );
    

    ##### get Seattle Seq data if in INFO field
    #FG:functionGVS FD:functionDBSNP GM:accession GL:geneList AAC:aminoAcids PP:proteinPosition 
    #CDP:cDNAPosition PP:polyPhen CP:scorePhastCons CG:consScoreGERP CADD:scoreCADD 
    #AA:chimpAllele RM:repeatMasker RT:tandemRepeat CA:clinicalAssociation DSP:distanceToSplice 
    #KP:keggPathway CPG:cpgIslands tfbs:transFactorBindingSites PPI:ProteinProteinInteraction 
    #PAC:proteinAccession GS:granthamScore MR:microRNAs
    push(@out, "==");
    my @ssFields = qw(
      GL FG  FD  GM  AAC  PP  CDP  AA  CA  DSP  KP  TFBS  PPI  PAC  CADD  PH  CP  CG  GS  MR
    );
    for my $item (@ssFields) {
      my $value =  (exists $info_hash{$item}) ? $info_hash{$item} : $MISSING;
      $value =~ s/^V\$// if ($item eq "TFBS");
      push(@out, $value);
    }

    #####  Annovar
    my $gene_detail = $info_hash{"AAChange.refGene"};
    $gene_detail = $info_hash{"GeneDetail.refGene"} if (not $gene_detail or $gene_detail eq ".");

    $info_hash{"Func.refGene"} =~ s/\\x3b/;/g;
    push(@out, "==",
        $info_hash{"Gene.refGene"} || $MISSING,
        $info_hash{"Func.refGene"} || $MISSING,
        $info_hash{"ExonicFunc.refGene"} || $MISSING,
        $gene_detail || $MISSING,
    );
    my @fields = qw(
        SIFT_score  SIFT_pred
        Polyphen2_HDIV_score  Polyphen2_HDIV_pred  Polyphen2_HVAR_score  Polyphen2_HVAR_pred
        LRT_score  LRT_pred
        MutationTaster_score  MutationTaster_pred  MutationAssessor_score  MutationAssessor_pred
        FATHMM_score  FATHMM_pred
        PROVEAN_score  PROVEAN_pred
        VEST3_score
        CADD_raw  CADD_phred
        DANN_score
        fathmm-MKL_coding_score  fathmm-MKL_coding_pred
        MetaSVM_score  MetaSVM_pred  MetaLR_score  MetaLR_pred
        integrated_fitCons_score  integrated_confidence_value
        GERP++_RS
        phyloP7way_vertebrate  phyloP20way_mammalian
        phastCons7way_vertebrate  phastCons20way_mammalian
        SiPhy_29way_logOdds
    );
    for my $item (@fields) {
      my $value =  (exists $info_hash{$item}) ? $info_hash{$item} : $MISSING;
      push(@out, $value);
    }
    
    
    # snpEff puts all info in one key: "EFF".  Each transcript is separated by a comma:
    # EFF=DOWNSTREAM(MODIFIER||4279||790|PERM1|protein_coding|CODING|NM_001291366.1||1),
    # NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Gcc/Acc|A176T|611|PLEKHN1|protein_coding|CODING|NM_032129.2|6|1)
    #
    # And the portion within the parenthesis can have these values (actual possible values are in quotes):
    # ("MODIFIER/HIGH/MODERATE/LOW" | "MISSENSE/NONSENSE/SILENT" | seq_change | AA_change | ??? |
    #   Gene | "protein_coding" | "CODING/NON_CODING" | Transcript | ?? | ??)
    $info_hash{"EFF"} //= "";
    my @snpEff = split(",", $info_hash{"EFF"});
    my (@snpEffGenes, @snpEffTypes, @snpEffDetails);

    # go through each transcript
    for my $transcript (@snpEff) {
        my ($mutationType, $mutationDetails) = split(/\(/, $transcript);
        my @details = split(/\|/, $mutationDetails);
        push(@snpEffGenes, $details[5]) if $details[5];
        push(@snpEffTypes, $mutationType);
        push(@snpEffDetails, $transcript);
    }
    @snpEffGenes = ($MISSING) if not @snpEffGenes;

    push(@out, "==",
                  join(",", uniq(@snpEffGenes)),
                  join(",", uniq@snpEffTypes),
                  join(",", @snpEffDetails));
    
    
    
    ##### get population frequencies (from annovar)
    push(@out, "==");
    my @freqFields = qw(
      esp6500siv2_ea
      esp6500siv2_aa
      esp6500siv2_all
      1000g2015aug_all
      1000g2015aug_afr
      1000g2015aug_amr
      1000g2015aug_eas
      1000g2015aug_eur
      1000g2015aug_sas
      ExAC_ALL
      ExAC_AFR
      ExAC_AMR
      ExAC_EAS
      ExAC_FIN
      ExAC_NFE
      ExAC_OTH
      ExAC_SAS
      ExAC_nontcga_ALL
      ExAC_nontcga_AFR
      ExAC_nontcga_AMR
      ExAC_nontcga_EAS
      ExAC_nontcga_FIN
      ExAC_nontcga_NFE
      ExAC_nontcga_OTH
      ExAC_nontcga_SAS

    );
    for my $item (@freqFields) {
      my $value =  (exists $info_hash{$item}) ? $info_hash{$item} : $MISSING;
      push(@out, $value);
    }
    
    

    say join("\t", @out);
    
# end while loop
}

print STDERR "\n";

close FILE;



