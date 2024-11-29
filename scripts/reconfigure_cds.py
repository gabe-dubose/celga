#!/usr/bin/env python3

import pyfaidx


#load data
sequence_file = '../data/genome/GCF_000002985.6/cds_from_genomic.fna'
sequence_data = pyfaidx.Fasta(sequence_file, key_function = lambda x:x.split('[locus_tag=')[1].split(']')[0], duplicate_action="longest", read_long_names=True)
sequence_ids = list(sequence_data.keys())


for sequence_id in sequence_ids:
    sequence = sequence_data[sequence_id]
    print(f">{sequence_id}\n{sequence}")