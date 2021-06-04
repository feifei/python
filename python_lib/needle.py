#!/usr/bin/python

import sys
import os
import subprocess
from optparse import OptionParser

from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO


def main():
    parser = OptionParser()
    parser.add_option("-i", "--infile1", dest="infile",
                      help="First sequence file", metavar="FILE")
    parser.add_option("-i", "--infile2", dest="infile2",
                      help="Second sequence file", metavar="FILE")

    parser.add_option("-o", "--outfile_aln", dest="outfile_aln",
                      help="aligned sequences in emboss text format", metavar="FILE")


    (options, args) = parser.parse_args()

    infile1 = options.infile1
    infile2 = options.infile2
    outfile = options.outfile_aln
    if infile1 == None  or infile2 == None or outfile == None:
        print "use -h or --help to see the help"
        sys.exit(1)

    run_needle_file(infile1, infile2, outfile)


def run_needle_file(infile1, infile2, outfile):
    cline = NeedleCommandline(asequence= infile1, bsequence = infile2, gapopen = 10, gapextend = 0.5, outfile = outfile)
    child = subprocess.Popen(str(cline), 
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    child.stderr.close()
    child.wait()
    return outfile


def run_needle(in_record1, in_record2):
    pipe_a, pipe_b = "/tmp/seq_pipe_a", "/tmp/seq_pipe_b"
    os.mkfifo(pipe_a), os.mkfifo(pipe_b)
    try:
        cline = NeedleCommandline(asequence = pipe_a, bsequence = pipe_b, gapopen = 10, gapextend = 0.5, outfile="/dev/stdout")
        child = subprocess.Popen(str(cline), 
                                 stdin=subprocess.PIPE, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE,
                                 shell=True)
        with open(pipe_a, 'w') as handle_a:
            SeqIO.write(in_record1, handle_a, "fasta")
        with open(pipe_b, 'w') as handle_b:
            SeqIO.write(in_record2, handle_b, "fasta")
    
        # Close STDERR to prevent its buffer from filling up, since we're not reading it.
        child.stdin.close()
        child.stderr.close() 
        return child.stdout
    finally:
        os.unlink(pipe_a), os.unlink(pipe_b)


if __name__ == "__main__":
    main()
