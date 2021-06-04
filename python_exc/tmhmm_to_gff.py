#!/usr/bin/env python

''' Convert TMHMM results to gff3 output
'''
import os, re

import argparse
from collections import defaultdict
from handy import gene_to_scf_pos, aa_to_na_pos

parser = argparse.ArgumentParser(description='Convert tmhmm out to gff file')
parser.add_argument('tmhmm_out')
parser.add_argument('gff_file', help = "gff file which holds the gene position")
args = parser.parse_args()

tmhmm_out = args.tmhmm_out
gff_file = args.gff_file

basename, extension = os.path.splitext(tmhmm_out)

tmhmm_gff_file = basename + ".gff"

icolor = "250 250 210"
ocolor = "219 219 112"
hcolor = "93 71 139"


genes = defaultdict(list)
with open(gff_file, 'r') as gff_fh:
    for line in gff_fh:
        line = line.strip()
        if len(line) == 0:
            continue
        arr = line.split("\t")
        scfid = arr[0]
        start, end = map(int, arr[3:5])
        strand = arr[6]
        geneid = re.match("ID=(.*?);", arr[-1]).groups()[0]
        genes[geneid] = [scfid, start, end, strand]
        


flag = False
total = 1
icount = 0
ocount = 0
hcount = 0
with open(tmhmm_gff_file, 'w') as outh, open(tmhmm_out, 'r') as th:
    for line in th:
        line = line.strip()
        if len(line) == 0 or re.match("#", line):
            if flag:
                flag = False
                total += 1
                icount, ocount, hcount = 0, 0, 0                    
            continue
        
        flag = True
        geneid, source, type, begin, end = line.split()
        if type == "TMhelix":
            hcount += 1
            attribs = "id=%d.h%d;color=%s;" %(total, hcount, hcolor)
        elif type == "inside":
            icount += 1
            attribs = "id=%d.i%d;color=%s;" %(total, icount, icolor)
        elif type == "outside":
            ocount +=1
            attribs = "id=%d.o%d;color=%s;" %(total, ocount, ocolor)

        scfid, scf_begin, scf_end, strand = genes[geneid]
        
        na_begin, na_end = aa_to_na_pos(int(begin), int(end))
        begin, end = gene_to_scf_pos(strand, scf_begin, scf_end, na_begin, na_end)
        
        print >>outh, "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s" %(scfid, source, type, begin, end, strand, attribs)
    