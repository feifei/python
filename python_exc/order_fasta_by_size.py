#!/usr/bin/env python
'''
    Order fasta files from biggest to smallest
'''
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', help="Fasta file to extract seq from")
args = parser.parse_args()

fasta_file = args.infile
outfile = fasta_file + '.ordered'


seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

with open(outfile, 'w') as outh:
    for seqid in sorted(seq_dict, key=lambda k: len(seq_dict[k]), reverse=True):
        print("%s %d" %(seqid, len(seq_dict[seqid])))
        outh.write(">%s\n" %seqid)
        outh.write(str(seq_dict[seqid]))
        outh.write("\n")
        
