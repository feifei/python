#!/usr/bin/env python
''' Convert GFF3 file to GTF, 
    only tested on EupathDB output, which is not a standard GFF3 file
    
'''

import os, re
import argparse

parser = argparse.ArgumentParser(description='Convert gff3 file to gtf file')
parser.add_argument('gff_file')
args = parser.parse_args()

gff_file = args.gff_file

basename, extension = os.path.splitext(gff_file)

gtf_file = basename + ".gtf"

with open (gff_file, 'r') as inh, open(gtf_file, 'w') as outh:
    for line in inh:
        if line.startswith("#"):
            continue
        scfid, source, feature, start, end, score, strand, frame, attr = line.split("\t")
        if feature == "gene":
            gene_id = re.match("ID=(.*?);", attr).group(1)
            new_attr = "gene_id \"%s\";" %gene_id
            transcript_id = gene_id + "_t"
        elif feature == "exon":
            exon_id = re.match("ID=(.*?);", attr).group(1)
            exon_number = re.search("-(\d+)$", exon_id).group(1)
            new_attr = "gene_id \"%s\"; transcript_id \"%s\"; exon_number \"%s\"; exon_id \"%s\";" %(gene_id, transcript_id, exon_number, exon_id)
        else:
            continue
        
        
        print >>outh, "\t".join([scfid, source, feature, start, end, score, strand, frame, new_attr])
        if feature == "gene":
            new_attr = "gene_id \"%s\"; transcript_id \"%s\";" %(gene_id, transcript_id)
            print >>outh, "\t".join([scfid, source, "transcript", start, end, score, strand, frame, new_attr])
        
