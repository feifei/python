#!/usr/bin/env python
'''
    Given gene description, and write the subset sequences in fasta file
'''

import re, os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file', help = "Fasta to subset from")
parser.add_argument('desc', type = str, help = "Gene description to extract from")
args = parser.parse_args()

fasta_file = args.fasta_file
desc = args.desc

records_subset = []
with open(fasta_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        description = record.description
        if re.search(desc, description, re.I):
            records_subset.append(record)

basename, extension = os.path.splitext(fasta_file)
fasta_file_sub = basename + ".sub.fasta"

with open(fasta_file_sub, 'w') as outh:
    SeqIO.write(records_subset, outh, 'fasta')