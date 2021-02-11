#!/usr/bin/env python3

from Bio import SeqIO
import os
import re
import sys

# First input: File with names of genomes
# Second input: Name of folder that contains alignments.
# All alignments should end in .faa and there should be no other files in that folder
# The gene names in the alignments should follow the prodigal format: GENOME_[0-9]*

# Get list of genomes
genomes = {}
with open(sys.argv[1], 'r') as handle:
    for x in handle.readlines():
        genomes[x.strip()] = ''

# Read through alignments and fill in
for file in os.listdir(sys.argv[2]):
    prots = {}
    
    # Load in seqs
    with open(os.path.join(sys.argv[2], file), 'r') as handle:
        for fa in SeqIO.parse(handle, 'fasta'):
            id_new = re.sub('_[0-9]*$', '', fa.id)
            prots[id_new] = fa.seq
    
    prot_len = len(fa.seq)

    # Loop through genomes and fill in
    for key in genomes:
        if key in prots.keys():
            genomes[key] = genomes[key] + prots[key]
        else:
            genomes[key] = genomes[key] + '-'*prot_len

for key in genomes:
    print('>'+key)
    print(genomes[key])
