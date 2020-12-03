#!/usr/bin/env python

'''
    Input gff file, fasta file and convert to genbank format
    One scaffold per genbank
    Usage:
    gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
'''

# Didn't work for rRNA file, something to do with parent child construction migh not be correct
# import sys
# import os
#
# from Bio import SeqIO
# from Bio.Alphabet import generic_dna
# from BCBio import GFF
#
# def main(gff_file, fasta_file):
#     out_file = "%s.gb" % os.path.splitext(gff_file)[0]
#     fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
#     gff_iter = GFF.parse(gff_file, fasta_input)
#     SeqIO.write(gff_iter, out_file, "genbank")
#
# if __name__ == "__main__":
#     main(*sys.argv[1:])

