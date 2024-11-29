#!/usr/bin/env python3

from gnat import phylogenetics


alignment_dir = '../data/strictcut_alignments'
extension = 'fasta'
phylogeny_dir = '../data/strictcut_phylogenies'

phylogenetics.batch_phylogenies(alignment_dir, extension, phylogeny_dir)