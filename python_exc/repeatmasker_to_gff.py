#!/usr/bin/env python

import argparse
import os, re

""" Convert infernal RepeatMasker .out result file to gff3 format
    RepeatMasker -species Diplomonadida -parallel 10 -gff ../data/ssk.cns.fa
"""


parser = OptionParser("usage: %prog repeatmasker_out_file")

parser.add_option("-s", "--source",
                  dest = "source", default = "RepeatMasker",
                  help="Specify source of the results")
parser.add_option("-t", "--type",
                  dest = "type", default = "repeat",
                  help="Specify the type of the results")
parser.add_option("-c", "--color",
                  dest = "color", default = "131 139 131",
                  help="Specify what color records from this file should be diplayed. Either as --color 2 or --color 255 0 0")
parser.add_option("-f", "--fragSize",
                  dest = "fragsize", default = 4000000,
                  help="Specify fragment size used for runing RepeatMasker")
parser.add_option("-l", "--overlapLen",
                  dest = "overlaplen", default = 2000,
                  help="Specify the overlap len used for running RepeatMasker")

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("incorrect number of arguments")

result_file = args[0]
source = options.source
type = options.type
color = options.color
frag_size = options.fragsize
overlap_len = options.overlaplen

basename, extension = os.path.splitext(result_file)

gff_file = basename + ".gff3"



def pos_in_scaffold(pos_in_frag, frag_num):
    ''' Return the scaffold position from the fragment position '''
    
    pos_in_scaffold = pos_in_frag + (frag_size - overlap_len) * frag_num 
    return pos_in_scaffold 


num = 1
with open(gff_file, 'w') as gff_h:
    with open(result_file, 'r') as infile_h:            
        for line in infile_h:
            line = line.strip()
            if len(line) == 0 or not re.match("\d", line):
                continue

            arr = line.strip("*").split()
            sw_score, div_perc, del_perc, ins_perc, qry, qry_start, qry_end, _, strand, repeat, repeat_class, repeat_start, repeat_end, _, ID = arr
            
            if re.search("frag-", qry):
                qry, frag_num = qry.split("frag-")
                frag_num = int(frag_num)
            else:
                frag_num = 0
            
            qry_start, qry_end = int(qry_start), int(qry_end)
            length = abs(qry_end - qry_start) + 1
            qry_start = pos_in_scaffold(qry_start, frag_num)
            qry_end = pos_in_scaffold(qry_end, frag_num)
            
            
            ID = "%s.%d" % (ID, num)
            num += 1
            attributes = "ID=%s;repeat=%s;repeat_class=%s;div_perc=%s;del_per=%s;ins_perc=%s;len=%d;repeat_start=%s;repeat_end=%s;color=%s;" %(ID, repeat, repeat_class, div_perc, del_perc, ins_perc, length, repeat_start, repeat_end, color)
            
            type = repeat_class.split("/")[0]
            
            print >>gff_h, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t.\t%s" %(qry, source, type, qry_start, qry_end, sw_score, strand, attributes)
            

