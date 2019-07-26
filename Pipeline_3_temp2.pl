#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';
use Getopt::Long;  #process command line arguments
use List::MoreUtils qw(first_index);



#### Global Pipeline Variables
my $SOFTWARE  = "/media/Data/software";        
my $PIPELINE  = "/media/Data/pipeline/v3";     
my $PUBLIC_DB = "/media/Data/public_databases";

# we don't do any filtering at this stage anymore.  Filtering is done after the VCFs
# are created, so that we can change filtering parameters without rerunning the 
# entire pipeline, including submitting to Seattle Seq
#may also want to remove all HOMO REF calls using: "MIN(AC)>=1 && (QUAL>=50  || (QUAL>=30 && MAX(FMT/DP)>=8))"
#my $QUALITY_FILTER = "QUAL>=50  || (QUAL>=30 && MAX(FMT/DP)>=8)"; 
#my $QUALITY_FILTER = "QUAL>=50  ";
my $QUALITY_FILTER = "";
#my $PASS_FILTER = "--apply-filters '.,PASS,FILTER'";  # for targeted seq we allow FILTER
my $PASS_FILTER = "--apply-filters ''";  # dont filter for now until we better understand VQS Lod filters.  

# my $KIDNEY_GENES = "$PUBLIC_DB/GeneLists/KidneyGenes6.txt";
# my $MOUSE_GENES  = "$PUBLIC_DB/GeneLists/CAKUT_genes-mouse.txt";
# my $OMIM         = "$PUBLIC_DB/GeneLists/OMIM_genemap.txt";
# my $HI           = "$PUBLIC_DB/HI/HI_prediction.bed";
# my $HI_IMP       = "$PUBLIC_DB/HI/HI_prediction_with_imputation.bed";
# my $RVI          = "$PUBLIC_DB/RVI/Genic_Intolerance_Scores.txt";
my $HG19_FASTA   = "$PUBLIC_DB/1000G/human_g1k_v37.fasta.gz";
# my $GDI          = "$PUBLIC_DB/GDI/GDI_full_10282015.txt";
# my $EXAC         = "$PUBLIC_DB/Exac/exac_LOF_alleles_by_gene.txt";
# my $CONSTRAINT   = "$PUBLIC_DB/Exac/forweb_cleaned_exac_r03_march16_z_data_pLI.txt";

my $JAVAPATH = "./:/media/Data/software/jars/*:" .
  "/media/Data/software/jars/httpunit-1.7/lib/*:" .
  "/media/Data/software/jars/httpunit-1.7/jars/*:";
$JAVAPATH .= $PIPELINE;





#### Get and Validate the arguments
my ($vcf_file, $Population, $MAF_Filter, $bFinalAnnotationOnly);
my $bSubmitToSeaSeq = 0;

GetOptions (
    'file|f=s' => \$vcf_file,
    'pop|p:s' => \$Population,
    'maf|m:f' => \$MAF_Filter,
    'sea-seq|s' => \$bSubmitToSeaSeq,
    'final-only|o' => \$bFinalAnnotationOnly,
)  or &usage("");

# Validate
( -e $vcf_file) or &usage("File not found: $vcf_file");
$MAF_Filter //= 1; # assign if undefined, 0 is ok (Novel variants)
$Population ||= "T";
($Population =~ /^(e|a|t|ea|ae)$/i)  or &usage("Invalid population designation.");
($MAF_Filter >= 0 && $MAF_Filter <= 1) or &usage("MAF threshold must be between 0 and 1.");



#### Create output file names
# Make sure VCF file is properly "bgzip" zipped
my $raw_vcf_file = $vcf_file;

($vcf_file !~ /\.gz$/i) || ($vcf_file !~ /\.vcf$/i) or &usage("File must be .vcf or vcf.gz.");

# Zip VCF file if necessary
if ($vcf_file !~ /\.gz$/i) {
    say "BGZipping VCF file...";
    system "bgzip $vcf_file";
    $vcf_file = "$vcf_file.gz";
}
$raw_vcf_file  =~ s/\.vcf\.gz$//i;  # strip .vcf.gz extension

my $rare = "MAF$MAF_Filter";
$rare = "Novel" if ($MAF_Filter == 0);
$rare = "NoMAF"  if ($MAF_Filter == 1);

my $ss = $bSubmitToSeaSeq ? "SS_" : "";

my $filt_file  =     "$raw_vcf_file.temp.vcf";
my $vcftools_filt =  "$raw_vcf_file.temp2.vcf.gz";

my $ss_file  =       "$raw_vcf_file.SS.vcf.gz";
my $snpeff_outfile = "$raw_vcf_file.$ss" . "Ann_Eff.vcf.gz";
my $table_file =     "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . ".tsv";
my $function_file =  "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . "_func.tsv";
my $disorder_file =  "$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . "_func_disorders.tsv";
my $comphet_summary ="$raw_vcf_file.$ss" . "Ann_Eff_" . $rare . "_CompHets.tsv";
my $input_file =     $filt_file;

if ($bFinalAnnotationOnly) {
    $snpeff_outfile = $vcf_file;
    $table_file  =     "$raw_vcf_file" . "_$rare" . ".tsv";   
    $function_file =   "$raw_vcf_file" . "_$rare" . "_func.tsv";
    $disorder_file =   "$raw_vcf_file" . "_$rare" . "_func_disorders.tsv";
    $comphet_summary = "$raw_vcf_file" . "_$rare" . "_CompHets.tsv";
    goto FINAL_ANN;
}



say "\n== Starting Pipeline ==";
say "== Input VCF must be aligned to hg19/37.";
say "== VCF file: $vcf_file...\n";



#### Index, QC and Normalize the VCF file with bcftools and vcftools
# Index the file
my $time = &getTime;
if (not -e "$vcf_file.csi") {
    say "== $time: Indexing VCF file...";
    system "bcftools index -f $vcf_file";
    $? == 0 or die "Failed to index VCF file.  Ensure the VCF file has been zipped using bgzip (not gzip)";
}

# Quality control: filters support MIN, MAX, AVG.
$time = &getTime;
say "== $time: Cleaning and Normalizing (filtering no longer happens before VCF files are annotated)";

#open(FILE, "bcftools view --exclude-uncalled --include '$QUALITY_FILTER' --apply-filters '$PASS_FILTER' $vcf_file | ") or
#    die "Couldn't filter $vcf_file";
    
open(FILE, "bcftools view --exclude-uncalled --include '$QUALITY_FILTER' $PASS_FILTER $vcf_file | ") or
    die "Couldn't filter $vcf_file";

# Note: if normalizing fails, debug VCF file with vcf-validator from the vcftools software package
open(FILT, "| bcftools norm -m -any -f $HG19_FASTA >  $filt_file") or die "Normalizing failed.";

while (my $line = <FILE>) {
    print FILT &cleanVCFLine(\$line);
}

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



#### Annovar Annotations
# see comments at end of this file for Annovar update/download instructions

$time = &getTime;
say "== $time: Running Annovar... ";

# REMEMBER TO UPDATE THESE IN Pipeline_VCFtoTable.pl/Pipeline_ReorderColumns.pl IF YOU CHANGE THEM!! #
my $proto = "refGene,dbnsfp33a,clinvar_20180603," .
            "popfreq_max_20150413," .
            "avsnp147," .
            "gnomad_exome," .
            "gme";


#             "esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all," .
#             "1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas";

my $cmd = "perl $SOFTWARE/annovar/table_annovar.pl $input_file $SOFTWARE/annovar/humandb/ 
    -buildver hg19 
    -out $raw_vcf_file.annovar
    -protocol $proto
    -operation g,f,f,f,f,f,f
    -nastring . 
    -vcfinput 
    -remove 
    -verbose";
$cmd =~ s/\n//g;

system $cmd;
$? == 0 or die "table_annovar.pl script failed";



#### SnpEff Annotations
$time = &getTime;
say "\n== $time: Running snpEff (snpEff generates a chr error on nonstandard chr i.e. GL000209.1 which can be ignored)...";
system "java -jar $SOFTWARE/snpEff/snpEff.jar -c $SOFTWARE/snpEff/snpEff.config GRCh37.75 $raw_vcf_file.annovar.hg19_multianno.vcf | 
          bgzip > $snpeff_outfile";
$? == 0 or die "snpEff.jar failed";






####### Start the custom scripts for final annotation, including kidney genes, known variants, etc #####
## if the -o flag is set, start here
FINAL_ANN:

$snpeff_outfile =~ /\.gz$/ or die "Input file $snpeff_outfile must be compressed with bgzip";


#### Genotype Quality Control
$time = &getTime;

say "== $time: No filtering done because GQ and/or DP can't be assumed to be provided in VCF.";
# say "== $time: Filtering Individual genotypes by min Depth 4 (min GQ 10 is now removed " .
#     "because Baylor seq does not have GQ so all variants were being filtered out)...";
#system "vcftools --gzvcf $snpeff_outfile --recode --recode-INFO-all --minDP 4 --minGQ 10 -c |
# system "vcftools --gzvcf $snpeff_outfile --recode --recode-INFO-all --minDP 4 -c |
#           bcftools view --exclude-uncalled |
#           bgzip > $raw_vcf_file.temp.vcf.gz";

system "bcftools view --exclude-uncalled $snpeff_outfile |
          bgzip > $raw_vcf_file.temp.vcf.gz";
$? == 0 or die "Vcftools filtering failed.\n  Try runing vcftools without '-c' and/or vcf-validator to find out why.";



#### Pipeline_VCFtoTable.pl extract all the values out of the INFO column and 
## makes make them into columns in the output table.  Also perform MAF filtering
$time = &getTime;
say "== $time: Converting VCF file to Table...";
#system "perl $PIPELINE/Pipeline_VCFtoTable_temp2.pl -m $MAF_Filter -f $raw_vcf_file.temp.vcf.gz > $table_file";
system "perl $PIPELINE/Pipeline_VCFtoTable_temp2.pl -m $MAF_Filter -f $raw_vcf_file.temp.vcf.gz > $table_file";


$? == 0 or die "Pipeline_VCFtoTable.pl script failed";



# MAF filtering is now handled by Pipeline_VCFtoTable.pl based on Annovar's Max MAF
# which is less powerful since a population can't be specified
#say "== Creating annotation file of rare variants using MAF < $MAF_Filter...";
#system "perl $PIPELINE/Pipeline_SelectRare.pl  $table_file $Population $MAF_Filter > $rare_file";
#  $? == 0 or die "Pipeline_SelectRare.pl script failed";



### filter variants by Deleterious Functional Mutations
$time = &getTime;
say "== $time: Filtering by functional variants...";
system "perl $PIPELINE/Pipeline_SelectFunctional.pl  $table_file  >  $function_file";
  $? == 0 or die "Pipeline_SelectFunctional.pl script failed";
  


#### Create a summary of potential Compound Heterozygous variants
$time = &getTime;
say "== $time: Creating summary file of potential compound hets...";
system "python $PIPELINE/Pipeline_CompHet.py < $function_file > $comphet_summary";
  $? == 0 or die ""; 



#### Add our custom Gene Based annotations including Emerge, Exac, etc.
$time = &getTime;
say "== $time: Creating gene based annotations...";
system "perl $PIPELINE/Pipeline_AdditionalAnnotations.pl < $function_file > $disorder_file.temp";
  $? == 0 or die "";
  

#### reorder columns
$time = &getTime;
say "== $time: Reordering columns";
system "perl $PIPELINE/Pipeline_ReorderColumns.pl < $disorder_file.temp > $disorder_file";
  $? == 0 or die "Reordering columns failed.";


#### Split output file into smaller files so they can be opened in Excel
# creates 2 files 
# -known genes: to see if the gene is known by clinvar, omim, or our 
#  custom lists for kidney genes and eMerge genes
# -deleterious: filters out variants that have functional types
# that are less likely to be pathogenic (INTRAGENIC|SPLICE_SITE_REGION|ncRNA|non-coding)
$time = &getTime;
say "== $time: Splitting main annotation file into smaller files";
system "perl $PIPELINE/Pipeline_SplitFile.pl $disorder_file";
  $? == 0 or die "Splitting file failed.";

#### Cleanup
say "== $time: Removing temporary files";
unlink glob("$raw_vcf_file.annovar*");
unlink glob("$raw_vcf_file.temp*"); 

unlink "$disorder_file.temp";

unlink glob("$filt_file*");  
unlink glob("snpEff_*");
unlink glob("out.log");

# Delete file filtered by functional variants.
# Not needed because we immediately generate this file again, but with the 
# disorders columns added (kidney, mouse, etc.)
# however, I do remove the INFO column in the disorders file, to make it smaller,
# so if you need the info column, don't delete this file.
unlink  $function_file; 

# delete unfiltered VCF -> table file 
# I no longer delete this because sometimes we want to see a file that is not filtered
# by variant type.
# Tip: If you need to view this file, you can delete the INFO column in this file 
# to make it a lot smaller.
#unlink $table_file;   



$time = &getTime;
say "== $time: Finished Pipeline..";







################  Functions  ################

sub usage {
    say "\nError: \n$_[0]\n\n
      Usage: ./Pipeline.pl -m 0.01 -f VCF_file  \n\n
      MAF_threshold is between 0 and 1 (ex. 0.01 for frequency 1 in 100), default is 0.01 \n\n";
    exit 1;
}
# Choosing a population has been removed for now since we use Annovar's Pop_max value.
# removed from usage message:
#Population options are: \n  T (Default - All Populations), E (Europian), A (African) and EA (Europian & African)


sub getTime() {
  my $ellapsedSeconds = time() - $^T;
  my $min = int($ellapsedSeconds / 60);
  my $sec = $ellapsedSeconds % 60;
  return sprintf("%01dm%02ds", $min, $sec);
}



sub cleanVCFLine( $ ){
    # May, 2016: The version of GATK used by IGM uses alternate allele "<*:DEL>" (or just *) 
    # to indicate that a variant falls within an previous defined deletion. 
    # Since we normalize to one variant per line, this doesn't really give us any additional 
    # information.  So I remove these lines because they break Seattle Seq. 
    # More info on spanning deletions: https://www.broadinstitute.org/gatk/guide/article?id=6926

    my $INFO_COLUMN_INDEX = 7;  #Info column is the 8th column in a VCF file
    my $line = ${$_[0]};
    
    #Lack of a description field for filters breaks vcftools
    if ($line =~ /##FILTER=/ && $line !~ /,Description=/) {
        $line =~ s/>/,Description="">/;
    }
    
    return $line if ($line =~ /^#/);
    return "" if $line =~ /^hs37d5/;  #some vcf files specify hs37d5 as a CHR which breaks bcftools
  
    chomp $line;
    
    # Trailing tabs may break annotation engines
    #$line =~ s/\t$//;

    # - Change missing genotypes to dipoloid so Annovar doesn't give warnings  (.:.:. should be ./.:.:.)
    # (this is done automatically by bcftools now)
    #$line =~ s|\t\.:|\t\./\.:|g;


    my @fields = split("\t", $line);

    # - Replace spaces in INFO field because SeaSeq will convert them to tabs, messing up the VCF file
    $fields[$INFO_COLUMN_INDEX] =~ s/ /_/g;

    # This function was created when dealing with old, poorly formed Yale VCF files, 
    # but since now most of the VCF files we get are pretty clean, this is no
    # longer necessary.
    
    # - since a value is required, give a missing value if it's not there
    #$line =~ s/DESCRIPTION2;/DESCRIPTION2=\.;/g;
    #$line =~ s/AACONSERV100;/AACONSERV100=\.;/g;
    #$line =~ s/LISTSPECIES100;/LISTSPECIES100=\.;/g;
    #$line =~ s/NOSPECIES100;/NOSPECIES100=\.;/g;
    
    # - Yale data sometimes has a pseudo alternate allele that acts as a descriptor for an indel.
    #   It starts with ",I" or ",D" so look for this and remove this since it breaks bcftools normalization
    #$fields[4] =~ s/,[ID][^,\t]+//;

    # - If the alterante allele starts with I, replace I with the reference allele
    #$fields[4] =~ s/^I/$fields[3]/;

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

