#!/usr/bin/env python

'''
    Filter only 4mC and 6mA, and color them red and blue respectively
'''

import re, os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff_file')

args = parser.parse_args()
gff_file = args.gff_file
basename, extension = os.path.splitext(gff_file)
outfile = "%s.%s%s" %(basename, 'colored', extension)

with open(gff_file, 'r') as inh, open(outfile, 'w') as outh:
    for line in inh:
        if line.startswith("#"):
            continue
        line = line.strip()
        arr = line.split("\t")
        type = arr[2]
        if type == "m4C":
            col = 2 #red
        elif type == "m6A":
            col = 4 #blue
        else:
            continue
        
        print >>outh, line + ";color=%d" %col
        
