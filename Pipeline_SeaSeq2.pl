#!/usr/bin/perl -w

=head1 USAGE Seattle Seq Annotation
================================================================================

=head1 Submit a VCF file to Seattle Seq for annotation

 Required:
 -i --input [file]
 -o --output [file]



 Example:
 ============
 Pipeline_SeaSeq.pl --input VCFFile.vcf.gz --output VCFFile.vcf.gz


=cut

use strict;
use warnings;
use Getopt::Long;
use LWP;
use File::Fetch;

my ($IN_FILE, $OUT_FILE);

GetOptions (
    'input|i=s' => \$IN_FILE,
    'output|o=s' => \$OUT_FILE,
    #'skip|s' => \$SKIP_FILTERS,
);
-e $IN_FILE or die;

(my $TEMP_FILE = $IN_FILE) =~ s/\.gz$//;
$TEMP_FILE =~ s/\.vcf$//;
$TEMP_FILE = "$TEMP_FILE.ss.temp.vcf";

system "echo '# autoFile testAuto\n# compressAuto true' > $TEMP_FILE";
system "bcftools view $IN_FILE  >> $TEMP_FILE ";

# print "== Zipping file before submission to Seattle Seq.\n";
# system "gzip -f $TEMP_FILE";
# $TEMP_FILE = "$TEMP_FILE.gz";

print "   Uploading file: $TEMP_FILE... \n\n";

my $browser = LWP::UserAgent->new;
push @{$browser->requests_redirectable}, 'POST';
my $response = $browser->post(
  'http://snp.gs.washington.edu/SeattleSeqAnnotation138/BatchQueryServlet',
  [
    "GenotypeFile" => [$TEMP_FILE => $TEMP_FILE => "application/x-gzip"],
    "autoFile" => "testAuto",   #doesn't seem to do anything, still need to add "# autoFile testAuto" to vcf file
    "instance" => "first",
    "genotypeSource" => "FileInput",
    "autoFile" => "none",
    "directWriteFile" => "none",
    "geneData" => "NCBI",
    "EMail" => 'davidfasel@gmail.com',
    "referenceTag" => "reference",
    "genotypeTag" => "genotype",
    "outputFileFormat" => "original",
    "chromosomeColumn" => "1",
    "locationColumn" => "2",
    "referenceColumn" => "0",
    "allele1Column" => "3",
    "allele2Column" => "4",
    "indelOutputFileFormat" => "newerIndel",
    "fileFormat" => "VCFSNVsAndIndels",
    "outputFileFormatBoth" => "VCFOutBoth",
    "columns" => "sampleAlleles",
    "columns" => "allelesDBSNP",
    "columns" => "scorePhastCons",
    "columns" => "consScoreGERP",
    "columns" => "scoreCADD",
    "columns" => "chimpAllele",
    "columns" => "geneList",
    "columns" => "HapMapFreq",
    "HapMapFreqType" => "HapMapFreqMinor",
    "columns" => "hasGenotypes",
    "columns" => "dbSNPValidation",
    "columns" => "repeats",
    "columns" => "proteinSequence",
    "columns" => "cDNAPosition",
    "columns" => "polyPhen",
    "polyPhenType" => "polyPhenScore",
    "columns" => "clinicalAssociation",
    "columns" => "distanceToSplice",
    "columns" => "microRNAs",
    "columns" => "grantham",
    "columns" => "keggPathway",
    "columns" => "cpgIslands",
    "columns" => "tfbs",
    "columns" => "genomesESP",
    "columns" => "PPI",
    "gFetch.x" => "74",
    "gFetch.y" => "2",
    "gFetch" => "Display Genotypes"
  ],
  'Content_Type' => 'form-data'
);

#die "Error: ", $response->status_line, "\n" =>  unless $response->is_success;
if (not $response->is_success) {
    die "Error: ", $response->status_line, "\n";
}
if ( $response->content =~ /html/i) {
    print $response->content;
    die "Something went wrong, the response from Seattle Seq is not as expected.";
}

unlink $TEMP_FILE;

my $str = $response->content;
my ($dlurl, $purl) = split (",", $str);
print "Successfully uploaded file. \nProgressURL:\n$purl\nDownloadURL:\n$dlurl\n\n";

while (1) {
    $response = $browser->get($purl);
    my ($lines_complete, $lines_total) = split (",", ( split("\n", $response->content) )[0]);
    print "Lines completed: $lines_complete of $lines_total\n";

    last if ($lines_complete > 0 && $lines_complete == $lines_total);
    sleep 10;
}

sleep 5;
my $ff = File::Fetch->new(uri => "$dlurl");
### fetch the uri to cwd() ###
my $where = $ff->fetch() or die $ff->error;

rename $where, $OUT_FILE;



