#!/usr/bin/env python

''' translate nucleotide sequences to aa sequence in fasta format'''

import sys
import os
import subprocess
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    args = parser.parse_args()

    infile = args.infile
    basename, extension = os.path.splitext(infile)
    outfile = basename + ".aa" + extension
    records = []
    with open(infile, 'r') as inh, open(outfile, 'w') as outh:
        for record in SeqIO.parse(inh, 'fasta'):
            if len(record.seq) %3 != 0:
                print record.id, "not dividable by 3"
            record.seq = record.seq.translate(6)
            records.append(record)
        SeqIO.write(records, outh, 'fasta')


if __name__ == "__main__":
    main()
