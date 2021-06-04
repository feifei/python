#!/usr/bin/python

import sys
import subprocess
from optparse import OptionParser
# from re import *

# from Bio import SeqIO
# from Bio.Seq import Seq

def main():
	parser = OptionParser()
	parser.add_option("-i", "--infile", dest="infile",
	                  help="fasta format of sequences to be aligned", metavar="FILE")
	parser.add_option("-o", "--outfile_aln", dest="outfile_aln",
	                  help="fasta format of corresponding aligned sequences", metavar="FILE")
	
	
	(options, args) = parser.parse_args()

	infile = options.infile
	outfile = options.outfile_aln
	if infile == None or outfile == None:
		print "use -h or --help to see the help"
		sys.exit(1)

	run_kalign(infile, outfile)
	

def run_kalign(infile, outfile):
	cl = KalignCommandline(infile, outfile)
	child = subprocess.Popen(str(cl), 
	                 stdout=subprocess.PIPE,
	                 stderr=subprocess.PIPE,
	                 shell=True)
	
	# Close pipes to prevent them from filling up, since we're not reading from them.
	child.stdout.close()
	child.stderr.close()
	
	# wait for the subprocess to finish first
	child.wait()
	assert child.returncode == 0, "kalign failed"

class KalignCommandline:
	""" Commandline for geneconv.
	"""
	def __init__(self, infile, outfile):
		self.infile = infile
		self.outfile = outfile
		self.parameters = []

	def __str__(self):
		 return "kalign < %s >| %s -f fasta" % (self.infile, self.outfile)


if __name__ == "__main__":
    main()
