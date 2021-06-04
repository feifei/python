#!/usr/bin/python
''' Universal alignment format converter
'''

from optparse import OptionParser
from Bio import AlignIO

def main():
    parser = OptionParser("usage: %prog infile outfile")
    parser.add_option("-i", "--informat", dest="informat", default= "fasta",
                      help="Format of the input file")
    parser.add_option("-o", "--outformat", dest="outformat", default = "phylip",
                      help="Format of the output file")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("incorrect number of arguments")

    informat = options.informat
    outformat = options.outformat
    infile = args[0]
    outfile = args[1]
    
    format_converter(infile, outfile, informat, outformat)


def format_converter(infile, outfile, informat, outformat):
    inh = open(infile, 'r')
    outh = open(outfile, 'w')

    alignments = AlignIO.parse(inh, informat)            
    AlignIO.write(alignments, outh, outformat)
    
    inh.close()
    outh.close()

    
if __name__ == "__main__":
    main()
