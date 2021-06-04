#!/usr/bin/env python
'''
    Calculate gaps (Ns) in fasta file
'''

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys, re
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file', type=argparse.FileType('r'), default=sys.stdin)
args = parser.parse_args()


gaps = defaultdict(list)
chr_sizes = defaultdict(int)
for record in SeqIO.parse(args.fasta_file, 'fasta'):
    chrid = record.id.split("_")[0]
    seq = str(record.seq).upper()
    gaps[chrid] = [(m.start(), m.end()) for m in re.finditer('N+', seq)]
    chr_sizes[chrid] = len(seq)


tot = 0
gap_sizes = []
for chrid, l in sorted(gaps.items()):
    if not l:
        print(chrid, "gap free")
    else:
        for start, end in l:
            gap_size = (end - start + 1)
            gap_sizes.append(gap_size)
            print(chrid, start, gap_size)
        tot += len(l)

print("There are %d gaps in total, and total size of %d bp." %(tot, sum(gap_sizes)))
print(sorted(gap_sizes))

chr_sizes = chr_sizes.values()
print("Sorted sizes of chromosomes:")
print(sorted(chr_sizes, reverse=True))

chromosomes = [i for i in chr_sizes if i >100000]
print("%d chromosomes with a total size of %d bp" %(len(chromosomes), sum(chromosomes)))