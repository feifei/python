#!/usr/bin/python

import sys
import subprocess
from optparse import OptionParser

from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO


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

    run_muscle_file(infile, outfile)


def run_muscle_file(in_file, out_file, **kwargs):
    cline = MuscleCommandline(input= in_file, out=out_file, **kwargs)
    child = subprocess.Popen(str(cline), 
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"))
    child.stderr.close()
    child.wait()
    return out_file


def run_muscle(in_records, **kwargs):
    cline = MuscleCommandline(input="/dev/stdin", out="/dev/stdout", **kwargs)
    child = subprocess.Popen(str(cline), 
                             stdin=subprocess.PIPE, 
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE,
                             shell=True)
    SeqIO.write(in_records, child.stdin, "fasta")
    
    # Close STDERR to prevent its buffer from filling up, since we're not reading it.
    child.stdin.close()
    child.stderr.close() 

    # child.wait() # if the output is in a file, it's needed
    # returncode is set by wait() or  poll()
    # assert child.returncode == 0, "Muscle failed"

    # return list(SeqIO.parse(child.stdout, "fasta"))
    return child.stdout
    

if __name__ == "__main__":
    main()
