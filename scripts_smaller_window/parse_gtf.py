#!/usr/bin/env python3

#A function to parse attributes field of GTF file
def parse_gtf_attributes(field, tag):
    for index, element in enumerate(field):
        if tag in element:
            return index

#A function to parase a GTF file and extract gene positions
def parse_gtf(file_path, outfile, cds_only=True):
    #initialize dictionary
    parsed_gtf = []

    #read file
    with open(file_path, 'r') as infile:
        lines = infile.readlines()
    
    #iterate through lines
    for line in lines:
        #pass comment lines
        if line[0] != '#':
            line = line.strip()
            #check gene entries only
            feature_type = line.split('\t')[2]
            #get information field 
            attributes = line.split('\t')[8]
            attributes_list = attributes.split(';')
            #check if entry is a gene
            if feature_type == 'gene':
                #check if only coding sequences should be included
                if cds_only == True:
                    if 'gene_biotype "protein_coding"' in attributes:
                        #get chromosome id
                        chromosome_id = line.split('\t')[0]

                        #get gene id
                        gene_id = attributes_list[parse_gtf_attributes(attributes_list, "gene_id")]
                        #clean up gene_id tag
                        gene_id = gene_id.lstrip('gene_id "').rstrip('"')
                        
                        #get position
                        start_position = int(line.split('\t')[3])
                        stop_position = int(line.split('\t')[4])

                        #add entry
                        parsed_gtf.append(f"{chromosome_id},{gene_id},{start_position},{stop_position}\n")

                #if operation is to be performed on all genes
                else:
                    #get chromosome id
                    chromosome_id = line.split('\t')[0]

                    #get gene id
                    gene_id = attributes_list[parse_gtf_attributes(attributes_list, "gene_id")]
                    #clean up gene_id tag
                    gene_id = gene_id.lstrip('gene_id "').rstrip('"')
                    
                    #get position
                    start_position = int(line.split('\t')[3])
                    stop_position = int(line.split('\t')[4])

                    #add entry
                    parsed_gtf.append(f"{chromosome_id},{gene_id},{start_position},{stop_position}\n")
        
    #write to file
    with open(outfile, 'a') as outfile:
        outfile.write(f"chromosome,gene,start,stop\n")
        for entry in parsed_gtf:
            outfile.write(entry)

#run
parse_gtf('../data/genome/GCF_000002985.6/genomic.gtf',
          '../data/gtf_parsed.csv')