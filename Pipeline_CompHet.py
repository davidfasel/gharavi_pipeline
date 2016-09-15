#!/usr/bin/python

import re, sys


print("\t".join(["Sample","Gene","Variants","Types"]))
chromosome_index = 0
gene_index = 0
types_index = 0
gene_dict = {}
index_of_parenthesis = 0


first_line = sys.stdin.readline().rstrip()

if not first_line.startswith("CHROM"):
    raise Exception("Invalid Header: Does not start with CHROM")
    
fields = first_line.split("\t")
sample_index = fields.index("Sample_ID(GT)")
gene_index = fields.index("GENES")
types_index = fields.index("MostDeleterious")

for line in sys.stdin:
    line = line.rstrip();
    
    # Make a list of samples, genes, and types, getting the fields
    # by using the indexes found earlier.        
    fields = line.split("\t")
    
    samples_list = fields[sample_index].split(",")
    all_samples = fields[sample_index]
    samples_separated_bycomma = all_samples.split(",")
    samples_separated_bycomma = filter(None, samples_separated_bycomma)
    
    all_genes = fields[gene_index].split(",")
    all_genes = filter(None, all_genes)
    tipes = fields[types_index]
    
    for sample in samples_separated_bycomma:
        
        # Each sample is checked to have its parenthesis removed
        # since only the type name is desired, whereas inside the parenthesis
        # it indicates whether if it's heterozygous or homozygous.
        if "(" in sample:
            index_of_parenthesis = sample.index("(")
            samples_cutoff = sample[:index_of_parenthesis]
        
        else:
            samples_cutoff = sample
        
        # Checks if the sample is put as a key in the dictionary.
        # If not, add it.
        if samples_cutoff not in gene_dict.keys():
            gene_dict [samples_cutoff] = {}
        
        # Goes through every gene and checks if the gene is already assigned
        # to the sample previously
        for gene in all_genes:
            
            # If the gene was not previously assigned, it assigns it to the sample
            # and sets up a counter for which it is incremented
            # everytime the gene appears for the sample
            # and sets the type.
            if gene not in gene_dict[samples_cutoff].keys():
                gene_dict[samples_cutoff][gene] = [1, tipes]
            
            # If the gene was previously assigned to the sample
            # it increments the gene counter for that specific sample by 1
            # and also adds the type.
            else: 
                gene_dict[samples_cutoff][gene][0] += 1
                gene_dict[samples_cutoff][gene].append(tipes)
                
# This nested loop prints each sample with its gene, the counter of the gene, and the type.
for sample in gene_dict.keys():

    for gene in gene_dict[sample].keys():
        
            if gene_dict[sample][gene][0]>1:
                
                print (sample+"\t"+gene+"\t"+str(gene_dict[sample][gene][0])+"\t"
                +str(",".join(gene_dict[sample][gene][1:])))