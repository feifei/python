#!/usr/bin/env python

''' Convert GFF3 file into separate embl files with one for each contig '''

from collections import defaultdict
import argparse
import os, re

from Bio import SeqIO
from Bio.Seq import Seq

from artemis_string import *
    
parser = argparse.ArgumentParser(description='Convert gff annotation file to separate EMBL format')
parser.add_argument('embl')
parser.add_argument('gff_file')
args = parser.parse_args()

run_dir = args.embl
gff_file = args.gff_file


def gff3_to_embl(arr):
    cds_prefix = "FT   CDS             "
    prefix = "FT                   "
    s = ""
    start = int(arr[3])
    end = int(arr[4])
    strand = arr[6]
    infos = arr[8]
    features_dict = defaultdict(str)
    for info in infos.split(";"):
        key, value = info.split("=")
        features_dict[key] = value
    geneid = features_dict["ID"]
    description = filter_string(features_dict["description"])
    if strand == "-":
        s += "%scomplement(%d..%d)\n" %(cds_prefix, start, end)
    else:
        s += "%s%d..%d\n" %(cds_prefix, start, end)
    s += "%s/geneid=\"%s\"\n" %(prefix, geneid)
    s += "%s/gene=\"%s [%s]\"\n" %(prefix, description, geneid)
    s += "%s/product=\"%s\"" %(prefix, description)
    return s


done = []
with open(gff_file, 'r') as gff_fh:
    for line in gff_fh:
        if len(line.strip()) == 0 or re.match("#", line):
            continue
        arr = line.strip().split("\t")
        scfid = arr[0]
        scfid = scfid.replace("/", "_")
        embl_record = gff3_to_embl(arr)
        if scfid not in done:
            # Make new file
            new_file = embl_out_dir + "/" + scfid + ".embl"
            with open(new_file, 'w') as fh:
                print >>fh, embl_record
            done.append(scfid)
        else :
            # continue to write
            with open(new_file, 'a') as fh:
                print >>fh, embl_record
                
                
                
