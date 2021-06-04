#!/usr/bin/env python3

'''
    Convert modification sites in gff to circlize highlights

'''


import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff_file')
parser.add_argument('outfile')

args = parser.parse_args()
gff_file = args.gff_file
outfile = args.outfile

with open(gff_file, 'r') as inh, open (outfile, 'w') as outh:
    for line in inh:
        line = line.strip()
        if line.startswith("#"):
            continue
        arr = line.split("\t")
        scfid = arr[0]
        type = arr[2]
        if not re.search("chr", scfid, re.I) or type == "modified_base":
            continue
        
        start, end = arr[3:5]
        strand = arr[6]
        if type == "m4C":
            col = "red"
        elif type == "m6A":
            col = "blue"
        
        outh.write("\t".join([scfid, start, end, col]) + '\n')
