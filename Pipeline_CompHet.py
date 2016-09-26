#!/usr/bin/python
import re, sys

gene_dict = {}

# print the header
print("\t".join(["Sample","Gene","Variants","Types"]))


first_line = sys.stdin.readline().rstrip()
if not first_line.startswith("CHROM"):
    raise Exception("Invalid Header: Does not start with CHROM")
    
# get the index of columns from the names of the columns in the headers
fields = first_line.split("\t")
sample_index = fields.index("Sample_ID(GT)")
gene_index = fields.index("GENES")
types_index = fields.index("MostDeleterious")

# go through the file and create a dictionary of sample IDs and corresponding genes
for line in sys.stdin:
    line = line.rstrip();
            
    fields = line.split("\t")
    
    samples_list = fields[sample_index].split(",")
    samples_list = filter(None, samples_list)
    
    all_genes = fields[gene_index].split(",")
    all_genes = filter(None, all_genes)
    
    func_type = fields[types_index]
    
    
    for sample in samples_list:
        
        # Each sample is checked to have its parenthesis removed
        # since only the type name is desired, whereas inside the parenthesis
        # it indicates whether if it's heterozygous or homozygous.
        if "(" in sample:
            index_of_parenthesis = sample.index("(")
            sample = sample[:index_of_parenthesis]
        
        # Checks if the sample a key in the dictionary.
        # If not, create an empty dictionary for that sample.
        if sample not in gene_dict.keys():
            gene_dict[sample] = {}
        
        # Go through every gene and check if the gene is already assigned
        # to the sample previously
        for gene in all_genes:
            
            # If the gene was not previously assigned, it adds the gene and the functional
            # type, and sets the count to 1 
            if gene not in gene_dict[sample].keys():
                gene_dict[sample][gene] = [1, func_type]
            
            # If the gene was previously seen for the sample,
            # increment the gene counter and add the type.
            else: 
                gene_dict[sample][gene][0] += 1
                gene_dict[sample][gene].append(func_type)
                
# This loop prints each sample with its gene, the count of the gene, and the type.
for sample in gene_dict.keys():

    for gene in gene_dict[sample].keys():
        if gene_dict[sample][gene][0] > 1:
            
            print (sample + "\t" + gene + "\t" + str(gene_dict[sample][gene][0])+ "\t" +
                   str(",".join(gene_dict[sample][gene][1:])))