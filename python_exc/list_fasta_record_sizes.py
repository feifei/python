#!/usr/bin/env python
'''
    List the sizes of the fasta records
'''

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', help="Fasta file to extract seq from")
args = parser.parse_args()

fasta_file = args.infile

with open(fasta_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id
        size = len(record.seq)
        print("%s\t%d" %(seqid, size))
        