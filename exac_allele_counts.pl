use warnings;
use strict;
use Array::Utils qw(unique  intersect);
use List::MoreUtils qw(first_index);
use feature "say";

my @header = qw(Gene  LOF.rare  Trunc.rare  Frame.rare  Splice.rare  Mis.rarePoly>0.9  
                      LOF  Trunc  Frame  Splice);

my @trunc = qw(
  stop_gained
  transcript_ablation
  start_lost
  stop_gained
  stop_lost
  initiator_codon_variant
);

# transcript_ablation
#   start_lost
#   stop_gained
#   stop_lost
#   initiator_codon_variant


my @splice = qw(
  splice_acceptor_variant
  splice_donor_variant
);
#splice_region_variant

my @frame = qw( frameshift_variant );

my @lof = (@trunc, @splice, @frame);
my %Genes;

while (<>) {
    chomp;
    my @fields = split(/\t/);
    
    next if (/^#/ || $fields[6] ne "PASS");
    
    my ($CSQ, $AC_Adj, $AN_Adj);
    my @info_fields = split(/;/, $fields[7]);
    
    
    for my $item (@info_fields) {
      my ($key, $value) = split("=", $item);
      if ($key eq "AC_Adj") {$AC_Adj = $value }
      elsif ($key eq "AN_Adj") {$AN_Adj = $value }
      elsif ($key eq "CSQ") {$CSQ = $value }
    }
    next if not ($AC_Adj && $CSQ);
    
    my $ref_allele = $fields[3];
    my @alt_alleles = split(/,/, $fields[4]);
    my @alt_allele_counts = split(/,/, $AC_Adj);
    
    # look through each transcript and look for LOF type mutations
    my %Alleles;
    my (@variant_genes);
    for my $trans (split(/,/, $CSQ)) {
        # get the gene names and mutation types for the transcript 
        my @transcript_fields = split(/\|/, $trans, -1);
        my $allele = $transcript_fields[0];
        my $gene_list = $transcript_fields[14];
        my $cons_list = $transcript_fields[4];
        my $confidence = $transcript_fields[48];
        my $polyphen  = $transcript_fields[25];
        next if (not $gene_list || not $cons_list);

        push (@variant_genes, $gene_list);
        my @consequences = split(/\&/, $cons_list);
        
        # correct for the funky way that CSQ encodes a deletion (ex. "-" for deletion)
        if ($allele eq "-") { $allele = substr($ref_allele, 0, 1); }
        my $allele_index = first_index {$_ eq $allele} @alt_alleles ;
        
        # if CSQ allele not found in ALT column, it's an Indel. Try normalizing it by adding first base from Ref allele
        if ($allele_index == -1) {
            $allele = substr($ref_allele, 0, 1) . $allele;
            $allele_index = first_index {$_ eq $allele} @alt_alleles ;  
            
            if ($allele_index == -1) {
                print STDERR "line: $. $ref_allele  $allele  $_";
                exit;
                
                next;
            }
        }
        $Alleles{$allele}{'index'} = $allele_index;
        
        if($confidence ne "LC" && intersect(@lof, @consequences)) {
            $Alleles{$allele}{'isLOF'} = 1;
            if(intersect(@trunc, @consequences))     { $Alleles{$allele}{'isStop'} = 1; }
            elsif(intersect(@frame, @consequences))  { $Alleles{$allele}{'isFrameshift'} = 1; }
            elsif(intersect(@splice, @consequences)) { $Alleles{$allele}{'isSplice'} = 1; }
        }
        
        if ($polyphen =~ /prob/) {  #probably deleterious
            $Alleles{$allele}{'isDelMissense'} = 1;
        }
    }
    
    
    # go through the list of Genes mentioned in all transcripts, and count the number of Alleles
    # associated with a lof variant
    @variant_genes = unique(@variant_genes);
    for my $g (@variant_genes) {
        for my $allele (keys %Alleles) {
            next if not $Alleles{$allele};
            my $i = $Alleles{$allele}{'index'};
    
            # if a LOF mutation was found, count it, and count only the worst subtype        
            if ($Alleles{$allele}{'isLOF'}) {
                $Genes{$g}{'lof'} += $alt_allele_counts[$i];
        
                if    ($Alleles{$allele}{'isStop'})       { $Genes{$g}{'trunc'}      += $alt_allele_counts[$i]; }
                elsif ($Alleles{$allele}{'isFrameshift'}) { $Genes{$g}{'frameshift'} += $alt_allele_counts[$i]; }
                elsif ($Alleles{$allele}{'isSplice'})     { $Genes{$g}{'splice'}     += $alt_allele_counts[$i]; }
            }
            
            # repeat the above but for rare alleles (MAF < 0.01)  
            next if ($alt_allele_counts[$i] / $AN_Adj > 0.01);
            if ($Alleles{$allele}{'isLOF'}) {
                $Genes{$g}{'lof_rare'} += $alt_allele_counts[$i];
        
                if    ($Alleles{$allele}{'isStop'})       { $Genes{$g}{'trunc_rare'}      += $alt_allele_counts[$i]; }
                elsif ($Alleles{$allele}{'isFrameshift'}) { $Genes{$g}{'frameshift_rare'} += $alt_allele_counts[$i]; }
                elsif ($Alleles{$allele}{'isSplice'})     { $Genes{$g}{'splice_rare'}     += $alt_allele_counts[$i]; }
            }
            if ($Alleles{$allele}{'isDelMissense'}) {
                $Genes{$g}{'mis_rare'} += $alt_allele_counts[$i];
            }
        }
    }
    print STDERR "\rline: $." if $. % 1000 == 0;
    
    #last if $. > 300000;
}


say join("\t", @header);

for my $gene (sort keys %Genes) {
    my @out;
    push(@out, 
      $gene, 
      $Genes{$gene}{'lof_rare'} || 0,  
      $Genes{$gene}{'trunc_rare'} || 0, 
      $Genes{$gene}{'frameshift_rare'} || 0, 
      $Genes{$gene}{'splice_rare'} || 0,
      $Genes{$gene}{'mis_rare'} || 0,
      $Genes{$gene}{'lof'} || 0,  
      $Genes{$gene}{'trunc'} || 0, 
      $Genes{$gene}{'frameshift'} || 0, 
      $Genes{$gene}{'splice'} || 0,
    );   
    say join("\t", @out);
}



=USAGE

cd /media/Data/public_databases/Exac/
zcat ExAC.r0.3.sites.vep.normalized.vcf.gz | 
  perl /media/Data/david/pipeline/exac_allele_counts.pl > exac_LOF_alleles_by_gene.txt
  
More info about the ExAC file here:
less /media/Data/public_databases/Exac/README

