#!/usr/bin/env python

''' Convert SignalP results to gff3 output
'''
import os, re

import argparse
from collections import defaultdict
from handy import gene_to_scf_pos, aa_to_na_pos

parser = argparse.ArgumentParser(description='Convert signalp out to gff file')
parser.add_argument('signalp_out')
parser.add_argument('gff_file', help = "gff file which holds the gene position")
args = parser.parse_args()

signalp_out = args.signalp_out
gff_file = args.gff_file

basename, extension = os.path.splitext(signalp_out)

signalp_gff_file = basename + ".gff"

color = "139 131 134"
type = "signal_pep"

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
        


num = 0
with open(signalp_gff_file, 'w') as outh, open(signalp_out, 'r') as sh:
    for line in sh:
        line = line.strip()
        if line.startswith("# SignalP"):
            source = re.search("(SignalP-.*?) ", line).groups()[0]
        if re.search("SP='YES'", line):
            m = re.match("Name=(.*?)\s*SP='YES' Cleavage site between pos. (\d+) and (\d+): (.*?) D=(.*?) D-cutoff=.*?Networks=(.*?)$", line)
            geneid = m.group(1)
            scfid, scf_begin, scf_end, strand = genes[geneid]
            cl_pos = int(m.group(2))
            cl_pattern = m.group(4)
            score = m.group(5)
            network = m.group(6)
                
            num += 1
            na_begin, na_end = aa_to_na_pos(1, cl_pos) 
            begin, end = gene_to_scf_pos(strand, scf_begin, scf_end, na_begin, na_end)
            
            attribs = "id=c.%d;score=%s;cl_pattern=%s;network=%s;color=%s;" %(num, score, cl_pattern, network, color)
            print >>outh, "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s" %(scfid, source, type, begin, end, strand, attribs)

