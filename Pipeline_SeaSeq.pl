#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';

my $PIPELINE = "/media/Data/david/pipeline";
my $classpath = "./:/media/Data/software/jars/*:" .
  "/media/Data/software/jars/httpunit-1.7/lib/*:" .
  "/media/Data/software/jars/httpunit-1.7/jars/*:";
$classpath .= $PIPELINE;


my $vcf_file   = $ARGV[0];
my $out_file   = $ARGV[1];
my $raw_vcf_file = $vcf_file;
my $temp_file;

( -e $vcf_file) or die "VCF file $vcf_file not found";
$raw_vcf_file =~ s/\.gz$//;
$raw_vcf_file =~ s/\.vcf$//;
$temp_file = "$raw_vcf_file.SSInput.temp.vcf";

if ($vcf_file =~ /\.gz$/) {
    open(FILE, "gunzip -c $vcf_file |") || die "can't open pipe to $vcf_file";
}
else {
    open(FILE, "$vcf_file") or die "Couldn't open $vcf_file";
}
open(TEMP, ">$temp_file") or die "Couldn't open $temp_file";

# add a comment line to the file required by SeaSeq
print TEMP "# autoFile testAuto\n# compressAuto true\n";
while (my $line = <FILE>) {
    print TEMP $line;
}

# zip the temp file
say "== BGzipping file before submission to Seattle Seq.";
system "bgzip -f $temp_file";
$temp_file = "$temp_file.gz";

say "== Starting online submission to Seattle Seq.";
system "java -cp $classpath SubmitSeattleSeqAnnotationAutoJob $temp_file $out_file";
    $? == 0 or die;

unlink $temp_file;

close FILE;
close TEMP;
# does this location still exist after processed by ss? 866511