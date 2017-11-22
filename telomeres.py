#!/usr/bin/env python
'''
    Look for contigs with telemeric repeats 
'''
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
import re

parser = OptionParser(usage = "usage: %prog [options] genome_fasta_file")

parser.add_option("-r", "--repeat", dest="repeat", default="TAGGG",
                  help="Specify telemer pattern")
parser.add_option("-n", "--repeat_num", dest="repeat_num", default=4,  type="int",
                  help="Specify telemer repeat number to search")
                  
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")
    
genome_file = args[0]
telemer = options.repeat * options.repeat_num
telemer_rc = str(Seq(telemer).reverse_complement())

with open(genome_file, 'r') as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        scfid = record.id
        scf_seq = str(record.seq).upper()
        if telemer in scf_seq :
            a = [m.start() for m in re.finditer(telemer, scf_seq)]
            start = a[0]
            end = a[-1]
            num = len(a)
            print scfid, "+1", start, end, num
        if telemer_rc in scf_seq:
            a = [m.start() for m in re.finditer(telemer_rc, scf_seq)]
            start = a[0]
            end = a[-1]
            num = len(a)
            print scfid, "-1", start, end, num


