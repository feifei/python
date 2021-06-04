#!/usr/bin/env python3
'''
    Calculate ATGCN in fasta file
'''

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys, re
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file', type=argparse.FileType('r'), default=sys.stdin)
args = parser.parse_args()



tot = 0
As = 0
Ts = 0
Cs = 0
Gs = 0
Ns = 0
for record in SeqIO.parse(args.fasta_file, 'fasta'):
    chrid = record.id.split("_")[0]
    seq = str(record.seq).upper()
    As += seq.count("A")
    Ts += seq.count("T")
    Cs += seq.count("C")
    Gs += seq.count("G")
    Ns += seq.count("N")
    tot += len(seq)


print("A: %d" %As)
print("T: %d" %Ts)
print("C: %d" %Cs)
print("G: %d" %Gs)
print("N: %d" %Ns)
print("tot: %d" %tot)