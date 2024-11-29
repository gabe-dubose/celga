#!/usr/bin/env python3

import pyfaidx
import json
import os

#A function to check that there are not redundant ids in homologous groups
def check_homolgoous_groups(file):

    #initialize list to store grouped ids and error code
    grouped = []
    error = 0

    #load data
    with open(file, 'r') as infile:
        cluster_data = json.load(infile)

    #iterate through each cluster
    for cluster_id in cluster_data:
        sequence_ids = cluster_data[cluster_id]
        #iterate through each sequence id
        for sequence_id in sequence_ids:
            if sequence_id not in grouped:
                grouped.append(sequence_id)
            elif sequence_id in grouped:
                print(f"Error: sequence {sequence_id} in multiple groups")
                error +=1
    if error == 0:
        print('Homologous groups check: pass')

    
    

#A function to partition gene sequence clusters into individual files
# input[clusters]: A json file where each key is a cluster name and each value is a list of sequence ids in said cluster
# input[sequences]: A fasta file that contains the sequences to be split
#Note: function requires pyfaidx
#Note: function assumes that the header line for each sequence is in NCBI format
def partition_sequence_clusters(clusters, sequences, outdir):

    #load cluster file
    with open(clusters, 'r') as infile:
        cluster_data = json.load(infile)

    #load sequence file
    sequence_data = pyfaidx.Fasta(sequences)
    #store sequence ids
    sequence_ids = list(sequence_data.keys())

    #iterate through each cluster
    for cluster in cluster_data:
        
        #add file
        os.system(f"touch {outdir}/{cluster}.fasta")
        #open file
        with open(f"{outdir}/{cluster}.fasta", 'w') as outfile:
            genes = cluster_data[cluster]
            #iterate through each gene in the cluster
            for gene_id in genes:
                #get sequence
                gene_sequence = str(sequence_data[gene_id])
                #write to file
                outfile.write(f">{gene_id}\n{gene_sequence}\n")

#using cds that I reconfigured ids to match
sequences = '/media/gabe/hd01/projects/active/C_elegans_lifecycle_genetics/data/genome/GCF_000002985.6/reconfigured_cds.fasta'

#Strict cutoffs
print('Partitioning sequences for strict cutoff groups...')
homologous_groups_file_strict = '../data/homologous_groups_strict_cutoffs.json'
check_homolgoous_groups(homologous_groups_file_strict)
strict_outdir = '../data/homologous_groups_seqs_strict'
#partition clusters
partition_sequence_clusters(homologous_groups_file_strict, sequences, strict_outdir)


#Lenient cutoffs
print('Partitioning sequences for lenient cutoff groups...')
homologous_groups_file_lenient = '../data/homologous_groups_lenient_cutoffs.json'
check_homolgoous_groups(homologous_groups_file_lenient)
lenient_outdir = '../data/homologous_groups_seqs_lenient'
#partition clusters
partition_sequence_clusters(homologous_groups_file_lenient, sequences, lenient_outdir)