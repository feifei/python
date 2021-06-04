#!/usr/bin/python

import sys
import subprocess
from optparse import OptionParser

# Not really working yet....

def main():
    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="infile",
                      help="Phylip format of the alignment file", metavar="FILE")
    parser.add_option("-o", "--outfile", dest="outfile",
                      help="Nexus format of tree result", metavar="FILE")


    (options, args) = parser.parse_args()

    infile = options.infile
    outfile = options.outfile
    if infile == None or outfile == None:
        print "use -h or --help to see the help"
        sys.exit(1)

    run_raxmlHPC(infile, outfile)

def run_raxmlHPC(in_file, out_file, args):
    cmd = ['raxmlHPC', '-m', 'PROTCATLGF', 
            '-s', in_file, '-n', out_file, 
            '-f', 'a', '-x', '12345', '-p', '12345', '-N', '100'] \
            + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    cline.stdout.close()
    cline.wait()


def run_raxmlHPC_PTHREADS(in_file, out_file, args):
    cmd = ['raxmlHPC-PTHREADS', '-m', 'PROTCATLGF', 
            '-s', in_file, '-n', out_file, 
            '-f', 'a', '-x', '12345', '-p', '12345', '-N', '100', '-T', '8'] \
            + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    cline.stdout.close()
    cline.wait()



if __name__ == "__main__":
    main()
