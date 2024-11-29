#!/usr/bin/env python3

from gnat import misc_utils
from gnat import phylogenetics

#define paths
seqs_path = '../data/homologous_groups_seqs_lenient'
extension = 'fasta'
alignment_outdir = '../data/lenientcut_alignments'

phylogenetics.batch_alignments(seqs_path, extension, alignment_outdir)