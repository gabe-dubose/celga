#!/usr/bin/env python3

print('Loading libraaries')

import json
from collections import defaultdict

#### Notes ####

# Two criteria for inferring homologous groups:
#	1. The query sequence must be at least 30% similar across the whole lenght of the target sequence (query coverage >= 100) and an e-value of >= 1e-10.
#	2. The query sequence must be at least 20% similar across 70% of the target sequence (query coverage >=70%) and an e-value of >= 1e-5.

# Blast results format:
#	6 qseqid sseqid pident length evalue bitscore qlen slen

# Alternative splices of the same gene will be considered the same (combined) because the data I'm using can't reliably quantify alternative splices
###############

print('Loading functions')
#A function to parse sequennce headers from an NCBI fasta file to get locus tags
def load_locus_info(locus_info_file):
	#initialize dictionary to store info
	id_to_locus_conversions = {}
	
	#read file
	with open(locus_info_file, 'r') as infile:
		lines = infile.readlines()
		
	#iterate through lines
	for line in lines:
		line = line.strip().split(' ')
		#separate fields
		try:
			locus_tag = [line[i] for i in range(len(line)) if 'locus_tag' in line[i]][0].replace('[locus_tag=', '').replace(']', '')
			protein_id = [line[i] for i in range(len(line)) if 'protein_id' in line[i]][0].replace('[protein_id=', '',).replace(']', '')
			
			#add info to output dictionary
			id_to_locus_conversions[protein_id] = locus_tag

		#skip pseudogenes
		except:
			pass
			
	return id_to_locus_conversions

def combine_groups(groups):
    # Create a defaultdict to store groups with common elements
    combined_groups = defaultdict(list)

    # Iterate through each group
    for group_key, group_values in groups.items():
        # Check against existing combined groups
        found_in_combined = False
        for combined_key, combined_values in combined_groups.items():
            if any(val in combined_values for val in group_values):
                # If a common element is found, combine groups
                combined_groups[combined_key].extend(group_values)
                found_in_combined = True
                break

        if not found_in_combined:
            # If no common element is found, create a new combined group
            combined_groups[group_key] = group_values

    return combined_groups

#a function to iteratively perform the combine function until all merging is complete
def iterative_combine(groups_dictionary):
	
	#initialize variable to define if merging has converged
	converged = 0

	#initialize starting dictionary
	starting_dictionary = groups_dictionary

	#initialize counter
	counter = 0

	#while merging has not converged
	while converged == 0:

		#performe combine function
		new_combine = combine_groups(groups_dictionary)

		#check if new combine is different
		if new_combine != starting_dictionary:
			#if so, redefine starting dictionary to be new combine run at next iteration
			starting_dictionary = new_combine
			counter += 1
			print(f"Joining on iteration: {counter}")
		#if new dictionary is the same as previous one
		elif new_combine == starting_dictionary:
			converged = 1
			print(f"Joining converged after {counter} iterations.")

			#remove redundancy in groups
			for group in new_combine:
				new_combine[group] = list(set(new_combine[group]))

			return new_combine

#A function to load blast results and keep matches of met criteria
def get_homologous_groups(blast_results_file, percent_identity_cutoff, alignment_coverage_cutoff, evalue_cutoff, locus_info):
	
	#initialize dictionary to store hits
	homologous_groups = {}
	
	#initialize dictionary for file writing and reformatting
	homologous_groups_json = {}
	
	#load blast results
	with open(blast_results_file, 'r') as infile:
		lines = infile.readlines()
		
	#iterate through lines
	for line in lines:
		try:
			#separate fields
			query_id = line.split('\t')[0]
			target_id = line.split('\t')[1]
			percent_identity = float(line.split('\t')[2])
			alignment_length = float(line.split('\t')[3])
			evalue = float(line.split('\t')[4])
			query_length = float(line.split('\t')[6])
			target_length = float(line.split('\t')[7].strip())
			
			#caculate alignment coverage
			alignment_coverage = alignment_length / target_length

			#get query and target locus tags
			query_locus = locus_info[query_id]
			target_locus = locus_info[target_id]
			
			#ignore hits that are one locus vs another
			if query_locus != target_locus:
				#check percent identity
				if percent_identity >= percent_identity_cutoff:
					#check alignment coverage
					if alignment_coverage >= alignment_coverage_cutoff:
						#check e value
						if evalue <= evalue_cutoff:
							#if all checks pass, add hit (locus tag) to dictionary
							#check if either entry is already in dictionary
							if query_locus in homologous_groups:
								homologous_groups[query_locus].append(target_locus)
							elif target_locus in homologous_groups:
								homologous_groups[target_locus].append(query_locus)
							else:
								homologous_groups[query_locus] = [target_locus]

		#skip message lines and pseudogenes
		except:
			pass
	
	#reformat homologous groups directory
	i = 0
	for group in homologous_groups:
		#assemble group
		group_id = f"group_{i}"
		full_group = homologous_groups[group]
		full_group.append(group)
		full_group = list(set(full_group))

		#add to dictionary
		homologous_groups_json[group_id] = full_group
		i+=1
	
	print('Combining groups...')
	homologous_groups_json = iterative_combine(homologous_groups_json)

	return homologous_groups_json

print('Starting calling of homologous groups.')

#define file paths
loci_info_file_path = '../data/loci_information.txt'
psi_blast_results_file_path = '../data/psiblast_results_Cele_vs_Cele.tsv'

print('Loading id to locus conversions...')
#load id to locus conversions - this to combine alternatively spliced genes
id_to_locus_conversions = load_locus_info(loci_info_file_path)

#strict cutoffs
#load blast data
print('Working on strict cutoffs...')
homologous_groups_strict_cutoffs = get_homologous_groups(psi_blast_results_file_path, 30, 1, 1e-10, id_to_locus_conversions)

print('Writing to file.')
#write to file
homologous_groups_strict_cutoffs_outfile = '../data/homologous_groups_strict_cutoffs.json'
homologous_groups_strict_cutoffs_json = json.dumps(homologous_groups_strict_cutoffs, indent=4)

with open(homologous_groups_strict_cutoffs_outfile, 'w') as outfile:
	outfile.write(homologous_groups_strict_cutoffs_json)

print(f"n groups with strict cutoffs: {len(homologous_groups_strict_cutoffs)}")

#lenient cutoffs
#load blast data
print('Working on lenient cutoffs...')
homologous_groups_lenient_cutoffs = get_homologous_groups(psi_blast_results_file_path, 20, 0.7, 1e-5, id_to_locus_conversions)

#write to file
print('Writing to file.')
homologous_groups_lenient_cutoffs_outfile = '../data/homologous_groups_lenient_cutoffs.json'
homologous_groups_lenient_cutoffs_json = json.dumps(homologous_groups_lenient_cutoffs, indent=4)

with open(homologous_groups_lenient_cutoffs_outfile, 'w') as outfile:
	outfile.write(homologous_groups_lenient_cutoffs_json)

print(f"n groups with lenient cutoffs: {len(homologous_groups_lenient_cutoffs)}")
