#!/usr/bin/env python
'''
    Given seqids, and write the subset sequences in fasta file
'''
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', help="Fasta file to extract seq from")
parser.add_argument('seqids', help="A list of geneids, with , to separate")
args = parser.parse_args()

fasta_file = args.infile
seqids = args.seqids.split(",")
with open(fasta_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id
        if seqid in seqids:
            print ">%s" %seqid
            print record.seq
        
