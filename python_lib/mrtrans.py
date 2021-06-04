#!/usr/bin/python

import sys
from optparse import OptionParser
from re import *

from Bio import SeqIO
from Bio.Seq import Seq

def main():
	parser = OptionParser()
	parser.add_option("-a", "--infile_aa_aln", dest="infile_aa_aln",
	                  help="fasta format of AA alignment sequences", metavar="FILE")
	parser.add_option("-n", "--infile_nt", dest="infile_nt",
	                  help="fasta format of corresponding NT sequences", metavar="FILE")
	parser.add_option("-o", "--outfile_nt_aln", dest="outfile_nt_aln",
	                  help="fasta format of corresponding NT sequences", metavar="FILE")
	
	(options, args) = parser.parse_args()

	aafile = options.infile_aa_aln
	ntfile = options.infile_nt
	outfile = options.outfile_nt_aln
	if ntfile == None or aafile == None or outfile == None:
		print "use -h or --help to see the help"
		sys.exit(1)
	
	coalign_files(ntfile, aafile, outfile)

def coalign_files(nt_infile, aa_infile, nt_outfile):
	nt_inhandle = open(nt_infile, 'r')
	aa_inhandle = open(aa_infile, 'r')
	ntout_handle = open(nt_outfile, 'w')

	nt_records = list(SeqIO.parse(nt_inhandle, "fasta"))
	aa_records = list(SeqIO.parse(aa_inhandle, "fasta"))

	coalign_records(nt_records, aa_records)

	SeqIO.write(nt_records, ntout_handle, "fasta")

	nt_inhandle.close()
	aa_inhandle.close()
	ntout_handle.close()

def order_nt_by_aa(nt_records, aa_records):
    ''' Order nt_records after aa_records 
        Muscle doesn' support stable anymore, has to reorder nt seqs manually
    '''
    new_nt_records = []
    for aa in aa_records:
        aa_id = aa.id
        for nt in nt_records:
            if nt.id == aa_id:
                new_nt_records.append(nt)
                break
    return new_nt_records

def coalign_records(nt_records, aa_records):
    assert len(nt_records) == len(aa_records), "AA file and NT file have different number of records."
    nt_records = order_nt_by_aa(nt_records, aa_records)
    for nt_record, aa_record in zip(nt_records, aa_records):
        assert aa_record.id == nt_record.id, "Name mismatch between AA record '"+aa_record.id+"' and NT record '"+nt_record.id+"'."
        nt_record.seq = coalign_seq(nt_record.seq, aa_record.seq)

def coalign_seq(nt_seq, aa_seq):
	codons = batch_gen(nt_seq, 3)
	out_seq = Seq('', nt_seq.alphabet)
	for codon in coalign(codons, aa_seq):
		out_seq += codon
	return out_seq

def coalign(codons, aas):
	for aa in aas:
		if aa == '-': yield '---'
		else: yield codons.next()
	#yield codons.next() # Stop codon

def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]

if __name__ == "__main__":
    main()
