# CPMG Annotation Pipeline


#### Running the pipeline

0. Make sure that you have all of the database files needed and change the paths...  Make sure any perl scripts are executable using `chmod +x <file>`.
1. Run the initial pipeline:  

        $pipeline/
        asdfasdfdsf

    Use the -o (final annotation *O*nly) to skip the SNPEff/Annovar annotation and just do the final annotations:

        $pipeline/

2. One var per line, CDW, etc etc etc.  Special attention to column numbers that have to be changed for 1 var/line. (or fix in file).

3. 



#### Description of Annotation Files

file1.pl
file2.pl
etc...

#### Notes
**SNPEff**
This version of the pipeline has removed the SNPEff annotation because that's expected to be done as part of the GATK alignment pipeline on Cavatica.

Todo: enable or disable the SNPEff annotation with a flag.

**Annovar**
comments about the version of the DBs used

**Omim**
Comments about downloading and using the latest OMIM files.

**FASTA File**
The fasta file that we used came in a tar: 
GRCh38_primary_assembly_plus_ebv_alt_decoy_hla.fasta.tar. The added decoys of  various DNA viruses help the alignment run faster but may not be that meaningful for our pipeline since we aren't doing FASTQ->BAM alignment.
I think we got it from Refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

We could consider switching to the fasta that Broad uses in Terra (although it may be basically the same one):
https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0/
