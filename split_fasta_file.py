#!/usr/bin/env python

import os
from optparse import OptionParser
from Bio import SeqIO



parser = OptionParser("usage: %prog fasta_file num_to_split ")

(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("incorrect number of arguments")

infile = args[0]
num = int(args[1])
file_name, file_ext = os.path.splitext(infile)

count = 0
records = []
file_num = 0
with open(infile, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        count += 1
        records.append(record)
        if count % num == 0:
            file_num += 1
            outfile = "%s.%d%s" %(file_name, file_num, file_ext)
            with open(outfile, 'w') as outh:
                SeqIO.write(records, outh, 'fasta')
            records = []
            count = 0

if records:
    file_num += 1
    outfile = "%s.%d%s" %(file_name, file_num, file_ext)
    with open(outfile, 'w') as outh:
        SeqIO.write(records, outh, 'fasta')
    