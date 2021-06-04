#!/usr/bin/env python
'''
    Given seqids with start-end, and write the subset sequences in fasta file
'''

import re, os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file', help = "Fasta to subset from")
parser.add_argument('info', type = str, help = "Geneid:start:end")
args = parser.parse_args()

fasta_file = args.fasta_file
info = args.info
wanted_seqid, start, end = info.split(":")
start = int(start)
end = int(end)

with open(fasta_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id.split("|")[0]
        if seqid == wanted_seqid:
            if end > start:            
                wanted_seq = record.seq[start - 1: end]
            else:
                wanted_seq = record.seq[end - 1 : start].reverse_complement()
            print ">%s" %info
            print wanted_seq
        
