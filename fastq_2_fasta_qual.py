#!/usr/bin/env python
'''
    Convert fastq file to fasta and qual
'''
import os
from optparse import OptionParser
from Bio import SeqIO

def main():
    parser = OptionParser("usage: %prog fastq_file")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("incorrect number of arguments")

    infile = args[0]
    file_name, file_ext = os.path.splitext(infile)
    
    fasta_outfile = file_name + ".fa"
    qual_outfile = file_name + ".qual"
    
    SeqIO.convert(infile, "fastq", fasta_outfile, "fasta")
    SeqIO.convert(infile, "fastq", qual_outfile, "qual")


if __name__ == "__main__":
    main()
