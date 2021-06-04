#!/usr/bin/env python
'''
    Print given seqs from fasta file,
    sub header to before |
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
        seqid = record.id.split("|")[0]
        if seqid in seqids:
            continue
        else:
            print(">%s" %seqid)
            print(record.seq)
        
