#!/usr/bin/env python
'''
    Given seqids, and write the subset sequences in fasta file
'''

from Bio import SeqIO
from optparse import OptionParser

parser = OptionParser("usage: %prog fasta_file seqids")

(options, args) = parser.parse_args()

if len(args) <= 1:
    parser.error("incorrect number of arguments")

fasta_file = args[0]
seqids = args[1:]

with open(fasta_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id.split("|")[0]
        if seqid in seqids:
            print ">%s" %seqid
            print record.seq
        
