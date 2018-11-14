#!/usr/bin/env python
'''
    Convert gff3 annotation to tab file, for import into R
'''

import argparse
import os, re
from collections import defaultdict
from artemis_string import *

parser = argparse.ArgumentParser(description='Convert gff3 annotation to tab file')
parser.add_argument('gff_file')
args = parser.parse_args()

gff_file = args.gff_file

basename, extension = os.path.splitext(gff_file)
tab_file = basename + ".tab"

with open(gff_file) as gff_fh, open(tab_file, 'w') as outh:
    print >>outh, "Geneid\tName\tDescription\tChromosome"
    for line in gff_fh:
        line = line.strip()
        if len(line) == 0:
            continue
        arr = line.split("\t")
        if len(arr) != 9:
            continue
        scfid, _, type, start, end, score, strand, _, infos = arr
        if not re.search("gene", type):
            continue
            
        geneid = re.search("ID=(.*?);", infos).groups()[0]
        if geneid.startswith("gene"):
            geneid = geneid.lstrip("gene:")
        m = re.search("description=(.*?);", infos)
        description = m.groups()[0] if m else ""
        if re.search("Source", description):
            description = description.split(" [Source")[0]
        name_m = re.search(";Name=(.*?);", infos)
        name = name_m.groups()[0] if name_m else ""
        description = filter_string(description)
        if description.endswith("mRNA."):
            description = re.match("(.*?)\(.*?\).*?, mRNA\.$", description).groups()[0]
        if not description:
            if re.search("biotype", infos):
                description = re.search("biotype=(.*?);", infos).groups()[0]
            else:
                description = ""

        print >>outh, "\t".join(map(str, [geneid, name, description, scfid]))