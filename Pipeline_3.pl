#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments
use List::MoreUtils qw(first_index);

# Global Pipeline Variables
my $SOFTWARE  = "./software";              #"/media/Data/software";
my $PIPELINE  = ".";                       #"/media/Data/pipeline/v3";
my $PUBLIC_DB = "./public_databases";      #"/media/Data/public_databases";


#may also want to remove all HOMO REF calls using: "MIN(AC)>=1 && (QUAL>=50  || (QUAL>=30 && MAX(FMT/DP)>=8))"
my $QUALITY_FILTER = "QUAL>=50  || (QUAL>=30 && MAX(FMT/DP)>=8)"; 
#my $QUALITY_FILTER = "QUAL>=50  ";
#my $PASS_FILTER = "--apply-filters '.,PASS,FILTER'";  # for targeted seq we allow FILTER
my $PASS_FILTER = "--apply-filters ''";  # dont filter for now until we better understand VQS Lod filters.  

my $KIDNEY_GENES = "$PUBLIC_DB/GeneLists/KidneyGenes6.txt";
my $MOUSE_GENES  = "$PUBLIC_DB/GeneLists/CAKUT_genes-mouse.txt";
my $OMIM         = "$PUBLIC_DB/GeneLists/OMIM_genemap.txt";
my $HI           = "$PUBLIC_DB/HI/HI_prediction.bed";
my $HI_IMP       = "$PUBLIC_DB/HI/HI_prediction_with_imputation.bed";
my $RVI          = "$PUBLIC_DB/RVI/Genic_Intolerance_Scores.txt";
my $HG19_FASTA   = "$PUBLIC_DB/1000G/human_g1k_v37.fasta.gz";
my $GDI          = "$PUBLIC_DB/GDI/GDI_full_10282015.txt";
my $EXAC         = "$PUBLIC_DB/Exac/exac_LOF_alleles_by_gene.txt";
my $CONSTRAINT   = "$PUBLIC_DB/Exac/forweb_cleaned_exac_r03_march16_z_data_pLI.txt";

my $JAVAPATH = "./:/media/Data/software/jars/*:" .
  "/media/Data/software/jars/httpunit-1.7/lib/*:" .
  "/media/Data/software/jars/httpunit-1.7/jars/*:";
$JAVAPATH .= $PIPELINE;
my $time = &getTime;

# Get and Validate the arguments
my ($vcf_file, $Population, $MAF_Filter, $bFinalAnnotationOnly);
my $bSubmitToSeaSeq = 0;



#### set options from paramaters
GetOptions (
    'file|f=s' => \$vcf_file,
    'pop|p:s' => \$Population,
    'maf|m:f' => \$MAF_Filter,
    'sea-seq|s' => \$bSubmitToSeaSeq,
    'final-only|o' => \$bFinalAnnotationOnly,
)  or &usage("");



#### Validate the Parameters
( -e $vcf_file) or &usage("File not found: $vcf_file");
$MAF_Filter //= 0.01; # assign if undefined, 0 is ok (Novel variants)
$Population ||= "T";
($Population =~ /^(e|a|t|ea|ae)$/i)  or &usage("Invalid population designation.");
($MAF_Filter >= 0 && $MAF_Filter < 1) or &usage("MAF threshold must be between 0 and 1.");


#### Mmake sure VCF file is properly "bgzip" zipped
my $raw_vcf_file = $vcf_file;

($vcf_file !~ /\.gz$/i) || ($vcf_file !~ /\.vcf$/i) or &usage("File must be .vcf or vcf.gz.");

# Zip VCF file if necessary
if ($vcf_file !~ /\.gz$/i) {
    say "BGZipping VCF file...";
    system "bgzip $vcf_file";
    $vcf_file = "$vcf_file.gz";
}
$raw_vcf_file  =~ s/\.vcf\.gz$//i;  # strip .vcf.gz extension


say "\n== Starting Pipeline for $vcf_file... \n";
say "\n== Please note that the input VCF must be aligned to hg19/37. \n";

# Output file names
my $rare = ($MAF_Filter == 0) ? "Novel" : "MAF$MAF_Filter";
my $ss = $bSubmitToSeaSeq ? "SS_" : "";
my $filt_file  =     "$raw_vcf_file.temp.vcf";
my $ss_file  =       "$raw_vcf_file.SS.vcf.gz";
my $annovar_out =    "$raw_vcf_file"."_annovar";
my $snpeff_infile =  "$raw_vcf_file"."_annovar.hg19_multianno.vcf";
my $snpeff_outfile = "$raw_vcf_file.$ss" . "Ann_Eff".".vcf.gz";
my $vcftools_filt =  "$raw_vcf_file.temp2.vcf.gz";
my $annovar_log =    "$raw_vcf_file.annovar.log";
my $table_file =     "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . ".tsv";
my $function_file =  "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . "_func.tsv";
my $disorder_file =  "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . "_func_disorders.tsv";
my $input_file =     $filt_file;

if ($bFinalAnnotationOnly) {
    $snpeff_outfile = $vcf_file;
    $table_file  =    "$raw_vcf_file" . "_$rare" . ".tsv";   
    $function_file =  "$raw_vcf_file" . "_$rare" . "_func.tsv";
    $disorder_file =  "$raw_vcf_file" . "_$rare" . "_func_disorders.tsv";
    goto FINAL_ANN;
}


#### Index, QC and Normalize the VCF file with bcftools and vcftools
# Index the file
if (not -e "$vcf_file.csi") {
    say "== Indexing VCF file...";
    system "bcftools index -f $vcf_file";
    $? == 0 or die "Failed to index VCF file.  Ensure the VCF file has been zipped using bgzip (not gzip)";
}

# Quality control: filters support MIN, MAX, AVG.
say "== Filtering by $QUALITY_FILTER...";

#open(FILE, "bcftools view --exclude-uncalled --include '$QUALITY_FILTER' --apply-filters '$PASS_FILTER' $vcf_file | ") or
#    die "Couldn't filter $vcf_file";
    
open(FILE, "bcftools view --exclude-uncalled --include '$QUALITY_FILTER' $PASS_FILTER $vcf_file | ") or
    die "Couldn't filter $vcf_file";



#### if normalizing fails, run VCF file through vcf-validator
open(FILT, "| bcftools norm -m -any -f $HG19_FASTA >  $filt_file") or die "Normalizing failed.";

while (my $line = <FILE>) {
    print FILT &cleanVCFLine(\$line);
}

#system "bgzip $filt_file";
#$filt_file = "$filt_file.gz";

close FILE;
close FILT;




#### Submit the file to Seattle Seq, if -s flag is set
if ($bSubmitToSeaSeq) {
    $input_file = $ss_file;
    
    say "== Starting online submission to Seattle Seq.";
    system "perl $PIPELINE/Pipeline_SeaSeq.pl $filt_file $ss_file";
    $? == 0 or die;

    #zipped file from SS won't index with bcftools, so we have to rezip with bgzip
    system "gunzip $ss_file";
    $ss_file =~ s/\.gz$//;
    
    # Seattle seq may convert lines with alt allele as "DEL" to lines that just have "null".  These lines will
    # break downline annotations, so remove these lines.
    system "grep -v '^null' $ss_file | bgzip > $ss_file.keep.gz"; 
    unlink $ss_file;
    rename "$ss_file.keep.gz", "$ss_file.gz";
    
    #system "bgzip $ss_file";
}
else {
    $input_file = $filt_file;
}





#### Start Annovar and SnpEff Annotations
# see comments at end of this file for download instructions
$time = &getTime;
say "== $time: Running Annovar... ";

######  REMEMBER TO UPDATE THESE IN Pipeline_VCFtoTable.pl IF YOU CHANGE THEM!!  #######
my $proto = "refGene,dbnsfp30a,clinvar_20160302,popfreq_max_20150413,exac03,exac03nontcga,avsnp147," .
            "esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all," .
            "1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas";

my $cmd = "perl $SOFTWARE/annovar/table_annovar.pl $input_file $SOFTWARE/annovar/humandb/ 
    -buildver hg19 
    -out $annovar_out 
    -protocol $proto
    -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f
    -nastring . 
    -vcfinput 
    -remove 
    -verbose";
$cmd =~ s/\n//g;
system $cmd;
$? == 0 or die "table_annovar.pl script failed";


#### run snpEff
$time = &getTime;
say "\n== $time: Running snpEff (snpEff generates a chr error on nonstandard chr i.e. GL000209.1 which can be ignored)...";
system "java -jar $SOFTWARE/snpEff/snpEff.jar -c $SOFTWARE/snpEff/snpEff.config GRCh37.75 $snpeff_infile | bgzip > $snpeff_outfile";
#system "java -jar $SOFTWARE/snpEff/snpEff.jar -c $SOFTWARE/snpEff/snpEff.config hg19 $snpeff_infile | bgzip > $snpeff_outfile";
$? == 0 or die "snpEff.jar failed";






### if the -o flag is set, start here
FINAL_ANN:

$time = &getTime;

$snpeff_outfile =~ /\.gz$/ or die "Input file $snpeff_outfile must be compressed with bgzip";

say "== $time: Filtering Individual genotypes by min Depth 4 and min GQ 10...";
system "vcftools --gzvcf $snpeff_outfile --recode --recode-INFO-all --minDP 4 --minGQ 10 -c | bcftools view --exclude-uncalled | bgzip > $vcftools_filt";
$? == 0 or die "Vcftools filtering failed.\n  Try runing vcftools without '-c' and/or vcf-validator to find out why.";



#### Our perl scripts start here


#### Basically extract all the values out of the INFO column and make them into columns in the 
#### output table
say "== $time: Converting VCF file to Table...";
system "perl $PIPELINE/Pipeline_VCFtoTable.pl -m $MAF_Filter -f $vcftools_filt > $table_file";
$? == 0 or die "Pipeline_VCFtoTable.pl script failed";


# MAF filtering is now handled by Pipeline_VCFtoTable.pl based on Annovar's Max MAF
# which is less powerful since a population can't be specified
# say "== Creating annotation file of rare variants using MAF < $MAF_Filter...";
# system "perl $PIPELINE/Pipeline_SelectRare.pl  $table_file $Population $MAF_Filter > $rare_file";
#   $? == 0 or die "Pipeline_SelectRare.pl script failed";



### filter variants by Deleterious Functional Mutations
say "== Filtering by functional variants...";
system "perl $PIPELINE/Pipeline_SelectFunctional.pl  $table_file  >  $function_file";
  $? == 0 or die "Pipeline_SelectFunctional.pl script failed";

say "== Creating gene annotations...";
system "perl $PIPELINE/Pipeline_AdditionalAnnotations.pl < $function_file > $disorder_file";
  $? == 0 or die "";
  
# reorder columns
say "== Reordering columns";
# desired column order
my @columns = qw(
  CHROM  POS  ID  avsnp147  REF  ALT  QUAL  FILTER  FORMAT  GENOTYPES  
  Sample_ID(GT)  Samples_With_Variant  AD_Pass  Missing 
  MAF_max  MostDeleterious  CADD  polyPhen  PhastCons  GERPConsScore  Clinvar
  
  GENES  RVI  RVI%  HI  HI_imp  HI%  HI%_imp GDI  GDI_Phred  GDI_Damage
    KidDisorder KidComment KidInheritence MouseGene MouseTerm OMIM_Disorder Emerge
    LOFRare.al  TruncRare.al  FrameRare.al  SpliceRare.al  MisRare.al.Poly>0.9
    LOF.al  Trunc.al  Frame.al  Splice.al 
    Exp_LOF.var  N_LOF.var Z_LOF  pLI
     
  ==SeattleSeq  geneList  functionGVS  functionDBSNP  accession  aminoAcids  proteinPosition
    cDNAPosition  chimpAllele  clinicalAssociation  distanceToSplice  keggPathway  tfbs
    PPI  proteinAccession  granthamScore  microRNAs
    
  ==Annovar  GeneAnn FuncAnn Exonic DetailsAnn SIFT_score SIFT_pred Polyphen2_HDIV_score 
    Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred 
    MutationTaster_score MutationTaster_pred MutationAssessor_score MutationAssessor_pred 
    FATHMM_score FATHMM_pred PROVEAN_score PROVEAN_pred VEST3_score CADD_raw CADD_phred 
    DANN_score fathmm-MKL_coding_score fathmm-MKL_coding_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred 
    integrated_fitCons_score integrated_confidence_value GERP++_RS phyloP7way_vertebrate phyloP20way_mammalian 
    phastCons7way_vertebrate phastCons20way_mammalian SiPhy_29way_logOdds
  
  ==snpEFF  GeneEFF  FuncEFF  DetailsEFF
  
  ==Freq All_ESP EUR_ESP AFR_ESP All_1KG Afr_1KG Amr_1KG Eas_1KG Eur_1KG Sas_1KG 
    ALL_Exac AFR_Exac AMR_Exac EAS_Exac FIN_Exac NFE_Exac OTH_Exac SAS_Exac  
    ALL_Exac_nontcga  AFR_Exac_nontcga  AMR_Exac_nontcga  EAS_Exac_nontcga  FIN_Exac_nontcga  NFE_Exac_nontcga  OTH_Exac_nontcga  SAS_Exac_nontcga
);
#SVM_PROBABILITY  SVM_POSTERIOR

rename $disorder_file, "$disorder_file.unordered.tsv";
open (IN, "$disorder_file.unordered.tsv");
open (OUT, "> $disorder_file");

my $line = <IN>;
chomp $line;
my @headers = split(/\t/, $line);
my @indexes;
for my $col (@columns) {
  my $index = first_index{$col eq $_} @headers;
  die "Can't reorder columns, couldn't find column header $col" if $index == -1;
  push(@indexes, $index);
}

say OUT join("\t", @columns);

while (<IN>) {
  chomp;
  my @fields = split("\t");
  my @out;
  push(@out, $fields[$_]) for (@indexes);  
  say OUT join("\t", @out);  
}



   
say "== Removing temporary files";
unlink "$disorder_file.unordered.tsv";
unlink glob("$annovar_out*multianno*"), "$annovar_out.avinput"; 
unlink glob("$filt_file*");  
unlink glob("snpEff_*");
unlink glob("out.log");
unlink $vcftools_filt;

unlink  $function_file; # output after functional filtering
unlink $table_file;   # output table before func filtering

$time = &getTime;
say "== $time: Finished Pipeline.";


###  Functions  ###

sub usage {
    say "Error: \n$_[0]\n\n
      Usage: ./Pipeline.pl -m 0.01 -f VCF_file  \n\n
      MAF_threshold is between 0 and 1 (ex. 0.01 for frequency 1 in 100), default is 0.01 \n\n";
    exit 1;
}
# removed from Usage message
#Population options are: \n  T (Default - All Populations), E (Europian), A (African) and EA (Europian & African)

sub getTime() {
  my $ellapsedSeconds = time() - $^T;
  my $min = int($ellapsedSeconds / 60);
  my $sec = $ellapsedSeconds % 60;
  return sprintf("%01dm%02ds", $min, $sec);
}

# This function was created when dealing with old, poorly formed Yale VCF files, 
# but since now most of the VCF files we get are pretty clean, I only replace spaces in 
# the info file because Seattle Seq chokes on them.

# May, 2016: The version of GATK used by IGM uses alternate allele "<*:DEL>" (or just *) 
# to indicate that a variant falls within an previous defined deletion. 
# Since we normalize to one variant per line, this doesn't really give us any additional 
# information.  So I remove these lines because they break Seattle Seq. 
# More info on spanning deletions: https://www.broadinstitute.org/gatk/guide/article?id=6926

sub cleanVCFLine( $ ){
    my $line = ${$_[0]};
    
    #Lack of a description field for filters breaks vcftools
    if ($line =~ /##FILTER=/ && $line !~ /,Description=/) {
        $line =~ s/>/,Description="">/;
    }
    
    return $line if ($line =~ /^#/);
    return "" if $line =~ /^hs37d5/;
  
    chomp $line;
    
    # Trailing tabs may break annotation engines
    #$line =~ s/\t$//;

    # - Change missing genotypes to dipoloid so Annovar doesn't give warnings  (.:.:. should be ./.:.:.)
    #$line =~ s|\t\.:|\t\./\.:|g;

    # - since a value is required, give a missing value if it's not there
#     $line =~ s/DESCRIPTION2;/DESCRIPTION2=\.;/g;
#     $line =~ s/AACONSERV100;/AACONSERV100=\.;/g;
#     $line =~ s/LISTSPECIES100;/LISTSPECIES100=\.;/g;
#     $line =~ s/NOSPECIES100;/NOSPECIES100=\.;/g;

    my @fields = split("\t", $line);

    # - Yale data sometimes has a pseudo alternate allele that acts as a descriptor for an indel.
    #   It starts with ",I" or ",D" so look for this and remove this since it breaks bcftools normalization
    #$fields[4] =~ s/,[ID][^,\t]+//;

    # - If the alterante allele starts with I, replace I with the reference allele
    #$fields[4] =~ s/^I/$fields[3]/;

    # - Replace spaces in INFO field because SeaSeq will convert them to tabs, messing up the VCF file
    $fields[7] =~ s/ /_/g;
    
    

    return join("\t", @fields) . "\n";
}

=scratch area
# dbSNP142 downloaded from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/

# download/update databases for annovar
cd /media/Data/software/annovar
./annotate_variation.pl -downdb -buildver hg19 refGene humandb
./annotate_variation.pl -downdb -buildver hg19 -webfrom annovar clinvar_20160302 humandb
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_ea humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_aa humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/  


=cut

