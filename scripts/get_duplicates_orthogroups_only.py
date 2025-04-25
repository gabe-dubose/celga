#!/usr/bin/env python3

import pandas as pd
import json

# load orthogroup data
orthogroup_data = pd.read_csv('../data/OrthoFinder/Results_Apr24/Phylogenetic_Hierarchical_Orthogroups/N0.tsv', sep='\t')
# rows with NA indicate that orthologs are not identified in all species -> remove
orthogroup_data.dropna(inplace=True)

# exclude one-to-one-to-one orthologs
def count_genes(s):
    if pd.isna(s):
        return 0
    return len(str(s).split(','))

orthogroup_data = orthogroup_data[
    (orthogroup_data['cele'].apply(count_genes) > 1) |
    (orthogroup_data['dplex'].apply(count_genes) > 1) |
    (orthogroup_data['dmel'].apply(count_genes) > 1)
]

# load duplicate data
# dmel
with open('../data/dmel_homologous_groups_strict_cutoffs.json', 'r') as infile:
    dmel_homologous_groups = json.load(infile)

#cele
with open('../data/duplication/homologous_groups_strict_cutoffs.json', 'r') as infile:
    cele_homologous_groups = json.load(infile)

#dplex
with open('../data/dpl_psiblast_id-30_e-neg-10_cov-1_sequence_clusters.json', 'r') as infile:
    dplex_homologous_groups = json.load(infile)

# get id conversions
# dmel
def parse_gtf_attributes(field, tag):
    for index, element in enumerate(field):
        if tag in element:
            return index
            
def parse_dmel_gtf(gtf_file):
    
    #initialize dictionary to store info
    locus_to_protein_id_conversions = {}
    locus_to_gb_id_conversions = {}
    id_conversions = {}
    
    # load file
    with open(gtf_file, 'r') as infile:
        gtf = infile.readlines()

    for line in gtf:
        # pass header
        if line[0] != '#':
            line = line.strip()
            #check gene entries only
            feature_type = line.split('\t')[2]
            #get information field 
            attributes = line.split('\t')[8]
            attributes_list = attributes.split(';')

            # get locus
            locus = attributes_list[parse_gtf_attributes(attributes_list, "gene_id")]
            locus = locus.lstrip('gene_id "').rstrip('"')

            try:
                #get gene id
                gene_id = attributes_list[parse_gtf_attributes(attributes_list, "db_xref")]
                #clean up gene_id tag
                gene_id = gene_id.lstrip('db_xref "').rstrip('"').replace("FLYBASE:", "")
                if gene_id[0:4] == 'FBgn':
                    locus_to_gb_id_conversions[locus] = gene_id
            except:
                pass

            try:
                # get protein id
                prot_id = attributes_list[parse_gtf_attributes(attributes_list, "protein_id")]
                # clean up
                prot_id = prot_id.lstrip('protein_id "').rstrip('"')
                locus_to_protein_id_conversions[locus] = prot_id
            except:
                pass
                
    for locus in locus_to_gb_id_conversions:
        if locus in locus_to_protein_id_conversions:
            prot_id = locus_to_protein_id_conversions[locus]
            gb_id = locus_to_gb_id_conversions[locus]
            id_conversions[prot_id] = gb_id
    
    return id_conversions       

dmel_id_conversions = parse_dmel_gtf('../data/dmel_genome/GCA_000001215.4/genomic.gtf')

# get dmel homologous groups

dmel_homologous_groups = {}

for index, row in orthogroup_data.iterrows():
    # get genes
    dmel_genes = row['dmel'].split(',')
    # make sure there are duplicates
    if len(dmel_genes) > 1:
        # add group to homologous groups
        group = f"group_{index}"
        dmel_homologous_groups[group] = []
        # go through each gene and get id conversion, and add
        for gene in dmel_genes:
            # get gene id
            gene = gene.strip()
            # get conversion
            if gene in dmel_id_conversions:
                id = dmel_id_conversions[gene]
                #add to group
                dmel_homologous_groups[group].append(id)
# clean up
for group in list(dmel_homologous_groups.keys()):
    group_size = len(dmel_homologous_groups[group])
    if group_size < 4:
        del dmel_homologous_groups[group]

# write to file
with open("../data/dmel_homologous_groups_orthologs_only.json", "w") as json_file:
    json.dump(dmel_homologous_groups, json_file, indent=4)

# get dplex groups

dplex_homologous_groups = {}

for index, row in orthogroup_data.iterrows():
    # get genes
    dplex_genes = row['dplex'].split(',')
    # make sure there are duplicates
    if len(dplex_genes) > 1:
        # add group to homologous groups
        group = f"group_{index}"
        dplex_homologous_groups[group] = []
        # go through each gene and get id conversion, and add
        for gene in dplex_genes:
            # get gene id
            gene = gene.strip()
            dplex_homologous_groups[group].append(gene)

# clean up
for group in list(dplex_homologous_groups.keys()):
    group_size = len(dplex_homologous_groups[group])
    if group_size < 4:
        del dplex_homologous_groups[group]

# write to file
with open("../data/dplex_homologous_groups_orthologs_only.json", "w") as json_file:
    json.dump(dplex_homologous_groups, json_file, indent=4)

# get id conversions for cele
def parse_gtf_attributes(field, tag):
    for index, element in enumerate(field):
        if tag in element:
            return index

def parse_cele_gtf(gtf_file):

    # initialize dictionary to store conversions
    cele_id_conversions = {}
    
    # load file
    with open(gtf_file, 'r') as infile:
        gtf = infile.readlines()
        
    for line in gtf:
        # pass header
        if line[0] != '#':
            line = line.strip()
            #check gene entries only
            feature_type = line.split('\t')[2]
            #get information field 
            attributes = line.split('\t')[8]
            attributes_list = attributes.split(';')
            #check if entry is a gene
            # get gene id
            gene_id = attributes_list[parse_gtf_attributes(attributes_list, "gene_id")]
            gene_id = gene_id.lstrip('gene_id "').rstrip('"').replace("CELE_", '')
            # get genbank id
            for index, element in enumerate(attributes_list):
                if 'db_xref "GenBank:' in element:
                    genbank_id = attributes_list[index]
                    genbank_id = genbank_id.lstrip('db_xref "GenBank:').rstrip('"')
                    cele_id_conversions[genbank_id] = gene_id

    return cele_id_conversions
    
cele_id_conversions = parse_cele_gtf('../data/genome/GCF_000002985.6/genomic.gtf')

# get cele homologous groups
cele_homologous_groups = {}

for index, row in orthogroup_data.iterrows():
    # get genes
    cele_genes = row['cele'].split(',')
    # make sure there are duplicates
    if len(cele_genes) > 1:
        # add group to homologous groups
        group = f"group_{index}"
        cele_homologous_groups[group] = []
        # go through each gene and get id conversion, and add
        for gene in cele_genes:
            # get gene id
            gene = gene.strip()
            # get id
            id = cele_id_conversions[gene]
            cele_homologous_groups[group].append(id)

# clean up
for group in list(cele_homologous_groups.keys()):
    group_size = len(cele_homologous_groups[group])
    if group_size < 4:
        del cele_homologous_groups[group]

# write to file
with open("../data/cele_homologous_groups_orthologs_only.json", "w") as json_file:
    json.dump(cele_homologous_groups, json_file, indent=4)