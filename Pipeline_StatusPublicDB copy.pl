#!/usr/bin/perl
#use 5.014;
use warnings;
use strict;
use feature 'say';

my $query=$ARGV[0];
my $dbSNPFile=$ARGV[1];
my $esp6500File=$ARGV[2];
my $CHROMOSOME = $ARGV[3];
$CHROMOSOME = "chr$CHROMOSOME";
#my $File1K=$ARGV[3];

my (%dbsnp, %db1k, %esp6500);

#define what should be shown if field is missing (can be empty string if desired)
my $MISSING = ".";

# NOTE!  if any of these headers change, they need to be updated in Pipeline_SelectRare.pl
my @header = qw(
    CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT GENOTYPES  
    Samples_Having_Variant  Sample_ID(GT)  Depths
    ==dbSNP  138_ID  141_ID  141_REF  141_ALT  141_AF  AlleleStatus
    ==ESP6500  ID  REF  ALT  Tot_AF  EUR_AF  AFR_AF  AlleleStatus
    ==1000g2014oct  all  afr  amr  eas  eur  sas
);
#former 1k fields (but now we use annovar to search 1kG)
#1K_ID  1K_REF  1K_ALT  1K_TOT_AF  1K_ASN_AF  1K_AMR_AF  1K_AFR_AF  1K_EUR_AF  1k_AlleleSatatus  ==

# todo, make Seattle Seq optional?
my @seattleseq_header = qw(
    ==SeattleSeq
    functionGVS
    functionDBSNP
    accession
    geneList
    aminoAcids
    proteinPosition
    cDNAPosition
    polyPhen
    scorePhastCons
    consScoreGERP
    scoreCADD
    chimpAllele
    repeatMasker
    tandemRepeat
    clinicalAssociation
    distanceToSplice
    keggPathway
    cpgIslands
    tfbs
    PPI
    proteinAccession
    granthamScore
    microRNAs
);
push (@header, @seattleseq_header);
say join("\t", @header) if $CHROMOSOME eq "chr1";

#print "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINDO\tFORMAT\tSamples_Having_Variant\tSample_ID(GT)\tDepths\t==\tdbSNP141_ID\tdbSNP141_REF\tdbSNP141_ALT\tdbSNP141_AF\tdbSNP141_AlleleStatus\t==\t1K_ID\t1K_REF\t1K_ALT\t1K_TOT_AF\t1K_ASN_AF\t1K_AMR_AF\t1K_AFR_AF\t1K_EUR_AF\t1k_AlleleSatatus\t==\tESP6500_ID\tESP6500_REF\tESP6500_ALT\tESP6500_Tot_AF\tESP6500_EUR_AF\tESP6500_AFR_AF\tESP6500_AlleleStatus\n";
#\t==\tKidney_Disorder\tInheritence\tOMIM Disorder\t";

#"Gene\tLength\tTot_Splice_Sites\t==\t";
#print "\#Rare_Stop_Variants\t\#Rare_Stop_Minor_alleles\t#Rare_Splice_Variants\t\#Rare_Splice_Minor_alleles\t\#Rare_Frameshift_Variants\t\#Rare_Frameshift_Minor_alleles\t\#Rare_Misense_Del_Variants\t\#Rare_Missense_Del_Minor_alleles\t\#Rare_Missense_Neutral_Variants\t\#Rare_Missense_Unknown_Variants\t\#Rare_Synonymous_variants\t\#Rare_Other_variants\t==\t\#Rare_Def_Del_variants\t\#Rare_All_Del_variants\t\#Total_Rare_variants\t\#Common_Stop_variants\t\#Common_Splice_Variants\t\#Common_Frameshift_Variants\t\#Common_Missense_Del_Variants\t\#Common_Misssense_Neutral_Variants\t\#Commom_Missense_Unknown_Variants\t#Common_Synonymous_Variants\t\#Common_Other_Variants\t==\t\#Common_Def_Del_Variants\t\#Common_All_Del_Variants\t\#Common_Total_Variants\t\#All_Stop\t\#All_Splice\t\#All_Frameshift\t\#All_Missense_Del\t\#All_Missense_Neutral\t\#All_Missense_Unknown\t\#All_Synonymous\t#All_Others\t==\tTotal_Def_Del\tTotal_All_Del\tTotal_Variants\n";



#### Parses the dbSNP database into a hash ####################
open(FILE,"$dbSNPFile") or die "File not found: $dbSNPFile ";
say STDERR "   Loading dbSNP database into memory...";
while(my $line = <FILE>) {
    #skip comment lines
    next if $line =~ /^#/;
    my ($chr, $pos) = split(/\t/, $line);
    $chr = "chr$chr" if ($chr !~ m/chr/);
    #only add lines to the hash if they are for the current chromosome
    next if $CHROMOSOME ne $chr;
    # create the key (i.e. chr1_123456) and store the line in the hash
    my $chr_pos = $chr."_".$pos;
    $dbsnp{$chr_pos} = $line;
}
close FILE;


######################### This section parses the 1K Genome database into a hash ########################
# open(FP2,"$File1K") or die "File not found: $File1K ";
# while($line2=<FP2>)
# {
#     chomp($line2);
#     @bb=split(/\t/,$line2);
#     if($bb[0]!~m/chr/)
#     {
#         $chr_1k="chr"."$bb[0]";
#     }
#     else{$chr_1k="$bb[0]";}
#     $chrom2="$chr_1k"."_"."$bb[1]";
#     $db1k{$chrom2}=$line2;
# }
# close FP2;


#### Parses the ESP6500 database into a hash 
open(FILE, "$esp6500File") or die "File not found: $esp6500File";
say STDERR "   Loading ESP6500 database into memory...";
while(my $line = <FILE>) {
    next if $line =~ /^#/;  #skip comment lines
    my ($chr, $pos) = split(/\t/, $line);
    $chr = "chr$chr" if ($chr !~ m/chr/);  #append 'chr' to chromosome if it's not there
    #only add lines to the hash if they are for the current chromosome
    next if $CHROMOSOME ne $chr;
    # create the key (i.e. chr1_123456) and store the line in the hash
    my $chr_pos = $chr."_".$pos;
    $esp6500{$chr_pos} = $line;
}
close FILE;


#### Build the Output file
my @sample_ids;
open(FILE,"$query") or die "File not found: $query";
while(my $line = <FILE>)
{
    chomp($line);
    my @var=split(/\t/,$line);
    my ($chr, $pos) = @var;
    my ($chrpos, $dp_index, $ad_index);
    my %info_hash;

    # get the sample IDs
    
    if($line =~ m/^\#CHROM/){
        for(my $s = 9; $s < @var; $s++){
            push (@sample_ids, $var[$s]);
        }
    }

    #skip comment, blank, or bad lines
    next if ($line =~ /^#/ or $line =~ /^\s*$/ or $line =~ /^\t/);

    # convert chr and postion to chr_pos format (ex. chr1_123456)
    $chrpos = "chr$chr" if ($chr !~ /chr/);
    $chrpos = $chrpos."_".$pos;

    #convert Info field into HASH table
    my @info_field = split(/;/, $var[7]);
    for my $item (@info_field) {
      my ($key, $value) = split("=", $item);
      #special case where item has no value, in which case it's a flag
      $value //= "flag";
      $info_hash{$key} = $value;
    }

    # get the index of the depth (DP) value in the FORMAT field.  
    # Format field may not exist if this is just a reference of variants without genotypes (such as EXAC)
    my $format = ($var[8] or "");
    if ($format) {
        my @dp_col = split(/\:/,$var[8]);
        for(0 .. @dp_col-1) {
            $dp_index = $_ if ($dp_col[$_] eq "DP");
            $ad_index = $_ if ($dp_col[$_] eq "AD");
        }
    }

    # count the samples and get their ID's, Genotypes, and Depths
    my ($samp, $depth) = ("", "");
    my $num_samples = 0;
    my $genotypes = "";
    for (my $c=9; $c < @var; $c++) {
        $genotypes .= "$var[$c];";
        my @gt_fields = split(/\:/, $var[$c]);
        my @gt = split("[/\|]", $gt_fields[0]);

        if( $gt[0] ne "." && $gt_fields[0] ne "0/0" ){
            $num_samples++;
        }
        if ($gt[0] ne ".") {
            $samp .= $sample_ids[$c-9] . "($gt_fields[0]),";
            
            # get the depths, prefer Allelic depth (AD) over just regular depth (DP)
            if ($ad_index && $gt_fields[$ad_index] ne ".") {
                $depth .= "$gt_fields[$ad_index];";
            }
            elsif ($dp_index && $gt_fields[$dp_index] ne ".") {
                $depth .= "$gt_fields[$dp_index];";
            }
            else {
                $depth .= ".;";
            }
        }
    }

    my @alleles_query = split(/\,/, $var[4]);
    
    my @out = (
      $var[0], $var[1], $var[2], $var[3], $var[4], $var[5], $var[6], $var[7], $format,
      $genotypes,
      $num_samples,
      $samp,
      $depth
    );

    # get dbsnp138 which was put in by Annovar
    my $value =  (exists $info_hash{'snp138'}) ? $info_hash{'snp138'} : $MISSING;
    push (@out, "==");
    push(@out, $value);

    # Get Info from DBSNP file.  Example format:
    # Y	10003	rs375039031	A	C	.	.	RS=375039031;RSPOS=10003;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x05000000000d000002000100;WGT=1;VC=SNV;CFL;ASP
    if(exists $dbsnp{$chrpos}) {
        my @fields = split(/\t/, $dbsnp{$chrpos});
        my @info_dbsnp = split(";", $fields[7]);
        my @alleles_dbSNP = split(",", $fields[4]);
        my $dbsnp_freq = $MISSING;

        for(my $i=0; $i < @info_dbsnp; $i++)
        {
            my @info_dbsnp_cols = split(/\=/, $info_dbsnp[$i]);
            if($info_dbsnp_cols[0] eq "CAF")
            {
                $info_dbsnp_cols[1]=~s/\[|\]//g;
                $dbsnp_freq=$info_dbsnp_cols[1];
            }
        }

        my $al1="DiffAllele";
        for(my $a1 = 0; $a1 < @alleles_query; $a1++) {
            for(my $a2=0; $a2<@alleles_dbSNP; $a2++) {
                if($alleles_query[$a1] eq $alleles_dbSNP[$a2]) {
                    $al1 = "SameAllele";
                    last;
                }
            }
        }
        push (@out, ($fields[2], $fields[3], $fields[4], $dbsnp_freq, $al1));
    }
    else {
        push (@out, ("Novel", $MISSING, $MISSING, $MISSING, $MISSING));
    }



    ##### get Info from 1000 genomes file (this is done by Annovar now, which is faster) #####
#     if(exists $db1k{$chrpos})
#     {
#         @var_1k=split(/\t/,$db1k{$chrpos});
#         @info_1k=split(/\;/,$var_1k[7]);
#         @alleles_1k=split(/\,/,$var_dbsnp[4]);
#         @maf_1k=split(/\=/,$info_1k[1]);
#
#         $tot_freq=$asn_freq=$afr_freq=$amr_freq=$eur_freq=0;
#         for($k1=0;$k1<@info_1k;$k1++)
#         {
#             @info_1k_cols=split(/\=/,$info_1k[$k1]);
#             if($info_1k_cols[0] eq "AF")
#             {
#                 $tot_freq=$info_1k_cols[1];
#             }
#             if($info_1k_cols[0] eq "ASN_AF")
#             {
#                 $asn_freq=$info_1k_cols[1];
#             }
#             if($info_1k_cols[0] eq "AMR_AF")
#             {
#                 $amr_freq=$info_1k_cols[1];
#             }
#             if($info_1k_cols[0] eq "AFR_AF")
#             {
#                 $afr_freq=$info_1k_cols[1];
#             }
#             if($info_1k_cols[0] eq "EUR_AF")
#             {
#                 $eur_freq=$info_1k_cols[1];
#             }
#         }
#
#         #output values are:
#         # 1K_ID  1K_REF  1K_ALT  1K_TOT_AF  1K_ASN_AF  1K_AMR_AF  1K_AFR_AF  1K_EUR_AF  1k_AlleleSatatus
#         print "\t==\t$var_1k[2]\t$var_1k[3]\t$var_1k[4]\t$tot_freq\t$asn_freq\t$amr_freq\t$afr_freq\t$eur_freq";
#
#         $al2="SameAllele";
#         for($a3=0;$a3<@alleles_query;$a3++)
#         {
#             $found2=0;
#             for($a4=0;$a4<@alleles_1k;$a4++)
#             {
#                 if($alleles_query[$a3] eq $alleles_1k[$a4])
#                 {
#                     $found2++;
#                 }
#             }
#             if($found2 == 0){
#                 $al2="DiffAllele Found";
#             }
#         }
#         print "\t$al2";
#     }
#     else
#     {
#         print "\t==\tNovel\t\t\t\t\t\t\t\t";
#     }


    ##### get Info from Exome Server Project file #####
    push (@out, "==");
    if(exists $esp6500{$chrpos})
    {
        my @var_esp = split(/\t/, $esp6500{$chrpos});
        my @info_esp = split(/\;/, $var_esp[7]);
        my @alleles_esp = split(/\,/, $var_esp[4]);
        my ($ea_maf, $aa_maf, $all_maf);

        for(my $c3=0; $c3 < @info_esp; $c3++){
            my @esp_info_cols = split(/\=/, $info_esp[$c3]);
            if($esp_info_cols[0] eq "MAF"){
                my @freq_esp = split(/\,/, $esp_info_cols[1]);
                $ea_maf = $freq_esp[0] / 100;
                $aa_maf = $freq_esp[1] / 100;
                $all_maf = $freq_esp[2] / 100;
            }
        }

        #compare alternate alleles to see if any of them match
        my $alleleStatus = "DiffAllele";
        for my $allele (@alleles_query) {
          for my $allele_esp (@alleles_esp) {
            $alleleStatus = "SameAllele" if ($allele eq $allele_esp);
          }
        }

        push (@out, ($var_esp[2], $var_esp[3], $var_esp[4], $all_maf, $ea_maf, $aa_maf, $alleleStatus));
    }
    else {
        push (@out, ("Novel", $MISSING, $MISSING, $MISSING, $MISSING, $MISSING, $MISSING));
    }


    ##### get 1000 genomes data from INFO field
    push(@out, "==");
    my @db1kGenFields = qw(
      1000g2014oct_all
      1000g2014oct_afr
      1000g2014oct_amr
      1000g2014oct_eas
      1000g2014oct_eur
      1000g2014oct_sas
    );
    for my $item (@db1kGenFields) {
      my $value =  (exists $info_hash{$item}) ? $info_hash{$item} : $MISSING;
      push(@out, $value);
    }



    ##### get Seattle Seq data if in INFO field
    # todo: make this optional?
    #FG:functionGVS FD:functionDBSNP GM:accession GL:geneList AAC:aminoAcids PP:proteinPosition 
    #CDP:cDNAPosition PP:polyPhen CP:scorePhastCons CG:consScoreGERP CADD:scoreCADD 
    #AA:chimpAllele RM:repeatMasker RT:tandemRepeat CA:clinicalAssociation DSP:distanceToSplice 
    #KP:keggPathway CPG:cpgIslands tfbs PPI PAC:proteinAccession GS:granthamScore MR:microRNAs
    push(@out, "==");
    my @ssFields = qw(
      FG  FD  GM  GL  AAC  PP  CDP  PH  CP  CG  CADD  AA  RM  RT  CA  DSP  KP  CPG  TFBS  PPI  PAC  GS  MR
    );
    for my $item (@ssFields) {
      my $value =  (exists $info_hash{$item}) ? $info_hash{$item} : $MISSING;
      push(@out, $value);
    }

    print join("\t", @out);



    #end the line
    print "\n";
# end while loop
}
say STDERR "   Done.";
close FILE;


