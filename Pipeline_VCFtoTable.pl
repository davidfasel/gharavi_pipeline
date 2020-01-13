#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments
use List::MoreUtils qw(uniq);
use List::Util 'max';
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

# NOTE!  if any of these headers change, they need to be updated in Pipeline_ReorderColumns.pl
my @header = qw(
    CHROM  POS  ID  avsnp150  REF  ALT  QUAL  FILTER  INFO  FORMAT GENOTYPES 
    Sample_ID(GT) Samples_With_Variant  Hets  Homs  AD_Pass  Missing  MAF_max  
        
    GENES

    ==snpEFF  snpEffGene  Annotation  Impact LOF
    
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
    G_AF 
    G_AF_raw 
    G_AF_male 
    G_AF_female 
    G_AF_afr 
    G_AF_ami 
    G_AF_amr 
    G_AF_asj 
    G_AF_eas 
    G_AF_fin 
    G_AF_nfe 
    G_AF_oth 
    G_AF_sas 
    GME_AF 
    GME_NWA 
    GME_NEA 
    GME_AP 
    GME_Israel 
    GME_SD 
    GME_TP 
    GME_CA 

    ==GeneScores
);
        
    # gnomad already includes 1000 Genomes and ESP
    # All_ESP  EUR_ESP  AFR_ESP     
    # All_1KG  Afr_1KG  Amr_1KG  Eas_1KG  Eur_1KG  Sas_1KG
         
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
    
    my %output;
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
    my @MAFs = (qw /   AF
                       AF_raw
                       AF_male
                       AF_female
                       AF_afr
                       AF_ami
                       AF_amr
                       AF_asj
                       AF_eas
                       AF_fin
                       AF_nfe
                       AF_oth
                       AF_sas
                       GME_AF
                       GME_NWA
                       GME_NEA
                       GME_AP
                       GME_Israel
                       GME_SD
                       GME_TP
                       GME_CA /);
                           
    my @maxMAF = (0);
    for my $maf (@MAFs) {
      if (looks_like_number $info_hash{$maf}) {
        push(@maxMAF, $info_hash{$maf});
      }
    }

    my $MAF_max = max(@maxMAF);
    
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
    
  
    
    $output{'CHROM'} = $var[0];
    $output{'POS'} = $var[1];
    $output{'ID'} = $var[2];
    $output{'avsnp150'} = $info_hash{"avsnp150"} || $MISSING;
    $output{'REF'} = $var[3];
    $output{'ALT'} = $var[4];
    $output{'QUAL'} = $var[5];
    $output{'FILTER'} = $var[6];
    $output{'INFO'} = $var[7];
    $output{'FORMAT'} = $format;
    $output{'GENOTYPES'} = $genotypes;
    $output{'Sample_ID(GT)'} = $samp;
    $output{'Samples_With_Variant'} = $num_samples;
    $output{'Hets'} = $hets;
    $output{'Homs'} = $homs;
    $output{'AD_Pass'} = $PASS_AD;
    $output{'Missing'} = $missing;
    $output{'MAF_max'} = $MAF_max;
    $output{'GENES'} = $gene;
    
    # $output{'SVM_PROBABILITY'} = &getItem(\%info_hash, 'SVM_PROBABILITY');
    # $output{'SVM_POSTERIOR'}   = &getItem(\%info_hash, 'SVM_POSTERIOR');
    
                
         


    #####  Annovar ########
    $output{'AAChange.refGene'} = &getItem(\%info_hash, 'AAChange.refGene');
    $output{'Func.refGene'} = &getItem(\%info_hash, 'Func.refGene');
    $output{'Gene.refGene'} = &getItem(\%info_hash, 'Gene.refGene');
    $output{'GeneDetail.refGene'} = &getItem(\%info_hash, 'GeneDetail.refGene');
    $output{'ExonicFunc.refGene'} = &getItem(\%info_hash, 'ExonicFunc.refGene');
    $output{'CADD_phred'} = &getItem(\%info_hash, 'CADD_phred');
    $output{'CADD_raw'} = &getItem(\%info_hash, 'CADD_raw');
    $output{'CADD_raw_rankscore'} = &getItem(\%info_hash, 'CADD_raw_rankscore');
    $output{'CLNALLELEID'} = &getItem(\%info_hash, 'CLNALLELEID');
    $output{'CLNDISDB'} = &getItem(\%info_hash, 'CLNDISDB');
    $output{'CLNDN'} = &getItem(\%info_hash, 'CLNDN');
    $output{'CLNREVSTAT'} = &getItem(\%info_hash, 'CLNREVSTAT');
    $output{'CLNSIG'} = &getItem(\%info_hash, 'CLNSIG');
    $output{'DANN_rankscore'} = &getItem(\%info_hash, 'DANN_rankscore');
    $output{'DANN_score'} = &getItem(\%info_hash, 'DANN_score');
    $output{'Eigen_coding_or_noncoding'} = &getItem(\%info_hash, 'Eigen_coding_or_noncoding');
    $output{'Eigen-PC-raw'} = &getItem(\%info_hash, 'Eigen-PC-raw');
    $output{'Eigen-raw'} = &getItem(\%info_hash, 'Eigen-raw');
    $output{'FATHMM_converted_rankscore'} = &getItem(\%info_hash, 'FATHMM_converted_rankscore');
    $output{'FATHMM_pred'} = &getItem(\%info_hash, 'FATHMM_pred');
    $output{'FATHMM_score'} = &getItem(\%info_hash, 'FATHMM_score');
    $output{'fathmm-MKL_coding_pred'} = &getItem(\%info_hash, 'fathmm-MKL_coding_pred');
    $output{'fathmm-MKL_coding_rankscore'} = &getItem(\%info_hash, 'fathmm-MKL_coding_rankscore');
    $output{'fathmm-MKL_coding_score'} = &getItem(\%info_hash, 'fathmm-MKL_coding_score');
    $output{'GenoCanyon_score'} = &getItem(\%info_hash, 'GenoCanyon_score');
    $output{'GenoCanyon_score_rankscore'} = &getItem(\%info_hash, 'GenoCanyon_score_rankscore');
    $output{'GERP++_RS'} = &getItem(\%info_hash, 'GERP++_RS');
    $output{'GERP++_RS_rankscore'} = &getItem(\%info_hash, 'GERP++_RS_rankscore');
    $output{'GTEx_V6p_gene'} = &getItem(\%info_hash, 'GTEx_V6p_gene');
    $output{'GTEx_V6p_tissue'} = &getItem(\%info_hash, 'GTEx_V6p_tissue');
    $output{'integrated_confidence_value'} = &getItem(\%info_hash, 'integrated_confidence_value');
    $output{'integrated_fitCons_score'} = &getItem(\%info_hash, 'integrated_fitCons_score');
    $output{'integrated_fitCons_score_rankscore'} = &getItem(\%info_hash, 'integrated_fitCons_score_rankscore');
    $output{'Interpro_domain'} = &getItem(\%info_hash, 'Interpro_domain');
    $output{'LRT_converted_rankscore'} = &getItem(\%info_hash, 'LRT_converted_rankscore');
    $output{'LRT_pred'} = &getItem(\%info_hash, 'LRT_pred');
    $output{'LRT_score'} = &getItem(\%info_hash, 'LRT_score');
    $output{'M-CAP_pred'} = &getItem(\%info_hash, 'M-CAP_pred');
    $output{'M-CAP_rankscore'} = &getItem(\%info_hash, 'M-CAP_rankscore');
    $output{'M-CAP_score'} = &getItem(\%info_hash, 'M-CAP_score');
    $output{'MetaLR_pred'} = &getItem(\%info_hash, 'MetaLR_pred');
    $output{'MetaLR_rankscore'} = &getItem(\%info_hash, 'MetaLR_rankscore');
    $output{'MetaLR_score'} = &getItem(\%info_hash, 'MetaLR_score');
    $output{'MetaSVM_pred'} = &getItem(\%info_hash, 'MetaSVM_pred');
    $output{'MetaSVM_rankscore'} = &getItem(\%info_hash, 'MetaSVM_rankscore');
    $output{'MetaSVM_score'} = &getItem(\%info_hash, 'MetaSVM_score');
    $output{'MutationAssessor_pred'} = &getItem(\%info_hash, 'MutationAssessor_pred');
    $output{'MutationAssessor_score'} = &getItem(\%info_hash, 'MutationAssessor_score');
    $output{'MutationAssessor_score_rankscore'} = &getItem(\%info_hash, 'MutationAssessor_score_rankscore');
    $output{'MutationTaster_converted_rankscore'} = &getItem(\%info_hash, 'MutationTaster_converted_rankscore');
    $output{'MutationTaster_pred'} = &getItem(\%info_hash, 'MutationTaster_pred');
    $output{'MutationTaster_score'} = &getItem(\%info_hash, 'MutationTaster_score');
    $output{'MutPred_rankscore'} = &getItem(\%info_hash, 'MutPred_rankscore');
    $output{'MutPred_score'} = &getItem(\%info_hash, 'MutPred_score');
    $output{'phastCons100way_vertebrate'} = &getItem(\%info_hash, 'phastCons100way_vertebrate');
    $output{'phastCons100way_vertebrate_rankscore'} = &getItem(\%info_hash, 'phastCons100way_vertebrate_rankscore');
    $output{'phastCons20way_mammalian'} = &getItem(\%info_hash, 'phastCons20way_mammalian');
    $output{'phastCons20way_mammalian_rankscore'} = &getItem(\%info_hash, 'phastCons20way_mammalian_rankscore');
    $output{'phyloP100way_vertebrate'} = &getItem(\%info_hash, 'phyloP100way_vertebrate');
    $output{'phyloP100way_vertebrate_rankscore'} = &getItem(\%info_hash, 'phyloP100way_vertebrate_rankscore');
    $output{'phyloP20way_mammalian'} = &getItem(\%info_hash, 'phyloP20way_mammalian');
    $output{'phyloP20way_mammalian_rankscore'} = &getItem(\%info_hash, 'phyloP20way_mammalian_rankscore');
    $output{'Polyphen2_HDIV_pred'} = &getItem(\%info_hash, 'Polyphen2_HDIV_pred');
    $output{'Polyphen2_HDIV_rankscore'} = &getItem(\%info_hash, 'Polyphen2_HDIV_rankscore');
    $output{'Polyphen2_HDIV_score'} = &getItem(\%info_hash, 'Polyphen2_HDIV_score');
    $output{'Polyphen2_HVAR_pred'} = &getItem(\%info_hash, 'Polyphen2_HVAR_pred');
    $output{'Polyphen2_HVAR_rankscore'} = &getItem(\%info_hash, 'Polyphen2_HVAR_rankscore');
    $output{'Polyphen2_HVAR_score'} = &getItem(\%info_hash, 'Polyphen2_HVAR_score');
    $output{'PROVEAN_converted_rankscore'} = &getItem(\%info_hash, 'PROVEAN_converted_rankscore');
    $output{'PROVEAN_pred'} = &getItem(\%info_hash, 'PROVEAN_pred');
    $output{'PROVEAN_score'} = &getItem(\%info_hash, 'PROVEAN_score');
    $output{'REVEL_rankscore'} = &getItem(\%info_hash, 'REVEL_rankscore');
    $output{'REVEL_score'} = &getItem(\%info_hash, 'REVEL_score');
    $output{'SIFT_converted_rankscore'} = &getItem(\%info_hash, 'SIFT_converted_rankscore');
    $output{'SIFT_pred'} = &getItem(\%info_hash, 'SIFT_pred');
    $output{'SIFT_score'} = &getItem(\%info_hash, 'SIFT_score');
    $output{'SiPhy_29way_logOdds'} = &getItem(\%info_hash, 'SiPhy_29way_logOdds');
    $output{'SiPhy_29way_logOdds_rankscore'} = &getItem(\%info_hash, 'SiPhy_29way_logOdds_rankscore');
    $output{'VEST3_rankscore'} = &getItem(\%info_hash, 'VEST3_rankscore');
    $output{'VEST3_score'} = &getItem(\%info_hash, 'VEST3_score');
    
    ##### get population frequencies (from annovar)
    $output{'G_AF'} = &getItem(\%info_hash, 'AF');
    $output{'G_AF_raw'} = &getItem(\%info_hash, 'AF_raw');
    $output{'G_AF_male'} = &getItem(\%info_hash, 'AF_male');
    $output{'G_AF_female'} = &getItem(\%info_hash, 'AF_female');
    $output{'G_AF_afr'} = &getItem(\%info_hash, 'AF_afr');
    $output{'G_AF_ami'} = &getItem(\%info_hash, 'AF_ami');
    $output{'G_AF_amr'} = &getItem(\%info_hash, 'AF_amr');
    $output{'G_AF_asj'} = &getItem(\%info_hash, 'AF_asj');
    $output{'G_AF_eas'} = &getItem(\%info_hash, 'AF_eas');
    $output{'G_AF_fin'} = &getItem(\%info_hash, 'AF_fin');
    $output{'G_AF_nfe'} = &getItem(\%info_hash, 'AF_nfe');
    $output{'G_AF_oth'} = &getItem(\%info_hash, 'AF_oth');
    $output{'G_AF_sas'} = &getItem(\%info_hash, 'AF_sas');
    $output{'GME_AF'}     = &getItem(\%info_hash, 'GME_AF');
    $output{'GME_NWA'}    = &getItem(\%info_hash, 'GME_NWA');
    $output{'GME_NEA'}    = &getItem(\%info_hash, 'GME_NEA');
    $output{'GME_AP'}     = &getItem(\%info_hash, 'GME_AP');
    $output{'GME_Israel'} = &getItem(\%info_hash, 'GME_Israel');
    $output{'GME_SD'}     = &getItem(\%info_hash, 'GME_SD');
    $output{'GME_TP'}     = &getItem(\%info_hash, 'GME_TP');
    $output{'GME_CA'}     = &getItem(\%info_hash, 'GME_CA');
    
    
    
    
    ####### snpEff ########
    
    # puts all info in three keys: "ANN", "LOF", NMD.  Each transcript is separated by a comma:

    ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
    #'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | 
    #Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | 
    #ERRORS / WARNINGS / INFO'">

    ##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | 
    #Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">

    ##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name 
    #| Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
 
    
    $info_hash{"ANN"} //= "";
    my @snpEff_ann = split(",", $info_hash{"ANN"});
    my (@snpEffGenes, @snpEffAnnotations, @snpEffAnnotation_Impact);

    # go through each transcript
    for my $transcript (@snpEff_ann) {
        my ($Allele, $Annotation, $Annotation_Impact, $Gene_Name, $Gene_ID, $Feature_Type, $Feature_ID, 
            $Transcript_BioType, $Rank, $HGVS_c, $HGVS_p, $cDNA_pos_length, $CDS_pos_length, $AApos_length, 
            $Distance, $ERRORS_WARNINGS_INFO) = split(/\|/,, $transcript);
        push(@snpEffGenes, $Gene_Name) if $Gene_Name;
        push(@snpEffAnnotations, $Annotation) if $Annotation;
        push(@snpEffAnnotation_Impact, $Annotation_Impact) if $Annotation_Impact;
    }
    @snpEffGenes = ($MISSING) if not @snpEffGenes;
    @snpEffAnnotations = ($MISSING) if not @snpEffAnnotations;
    @snpEffAnnotation_Impact = ($MISSING) if not @snpEffAnnotation_Impact;

    
    $output{'snpEffGene'} = join(",", uniq(@snpEffGenes));
    $output{'Annotation'} = join(",", uniq(@snpEffAnnotations));
    $output{'Impact'}     = join(",", uniq(@snpEffAnnotation_Impact));

    $output{'LOF'} = &getItem(\%info_hash, 'LOF');

     
    
    

    
    ##### output to file ########
    
    my @out;
    for my $h (@header) {
      #say STDERR $h, $output{$h};
      if ($h =~ /==/) {
          push(@out, "==");
      }
      else {
        #   if (not exists $output{$h}) {exit(1)};
          push(@out, $output{$h});
      }
    }
    

    say join("\t", @out);
    
# end while loop
}

print STDERR "\n";
close FILE;



sub getItem() {
  my %dict = %{$_[0]};
  my $item = $_[1];
#   if(not exists $dict{$item}) {print STDERR "missing $item;"; exit(1)}
  return (exists $dict{$item}) ? $dict{$item} : $MISSING;
}



