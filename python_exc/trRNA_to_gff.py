#!/usr/bin/env python

''' Make tRNA and rRNA gff file from
    tRNAscan-SE and rnammer output
'''

import argparse
import sys, os
import subprocess
import re
from Bio import SeqIO

def gff_record(id, scfid, source, feature, pos, score, strand, description, size):
    ''' Make a RNA gff3 record'''
        
    gene_attributes = "ID=%s;size=%d;description=%s" %(id, size, description)
    RNA_id = "rna_%s-1" %id
    RNA_attributes = "ID=%s;Parent=%s;size=%d;description=%s" %(RNA_id, id, size, description)
       
    if len(pos) > 2:
        start, intron_start, intron_end, end = pos
        exon1_id = "exon_%s-1" %id
        exon2_id = "exon_%s-2" %id
        exon1_attributes = "ID=%s;Parent=%s;description=exon" %(exon1_id, RNA_id)
        exon2_attributes = "ID=%s;Parent=%s;description=exon" %(exon2_id, RNA_id)
        record =  "\t".join(map(str, [scfid, source, "gene", start, end, score, strand, ".", gene_attributes])) + "\n"
        record += "\t".join(map(str, [scfid, source, feature, start, end, ".", strand, ".", RNA_attributes])) + "\n"
        record += "\t".join(map(str, [scfid, source, "exon", start, intron_start - 1, ".", strand, ".", exon1_attributes])) + "\n"
        record += "\t".join(map(str, [scfid, source, "exon", intron_end + 1, end, ".", strand, ".", exon2_attributes])) + "\n"
        
    else:
        start, end = pos
        exon_id = "exon_%s-1" %id
        exon_attributes = "ID=%s;Parent=%s;description=exon" %(exon_id, RNA_id)
        record =  "\t".join(map(str, [scfid, source, "gene", start, end, score, strand, ".", gene_attributes])) + "\n"
        record += "\t".join(map(str, [scfid, source, feature, start, end, ".", strand, ".", RNA_attributes])) + "\n"
        record += "\t".join(map(str, [scfid, source, "exon", start, end, ".", strand, ".", exon_attributes])) + "\n"
    
    return record


def process_RNA_file(RNA_file, source, prefix):
    global rno, tno
    
    records = ""
    with open(RNA_file, 'r') as th:
        for line in th:
            line = line.strip()
            if re.match("#", line) or len(line) == 0:
                continue
        
            if source == "rnammer":
                # gff2 
                scfid, _, feature, start, end, score, strand, _, description = line.split("\t")
            
                start, end = int(start), int(end)
                size = end - start + 1
                pos = [start, end]
                
                if description == "8s_rRNA":
                    description = "5S rRNA"                
                description = description.replace("s_", "S ")
                
                rno += 1
                id = "%s%04d" %(prefix, rno)
                        
            if source == "tRNAscan-SE":
                feature = "tRNA"
                if re.match("Sequence|Name|-", line) or len(line) == 0:
                    continue                
                    
                scfid, _, start, end, aa, anti_codon, intron_start, intron_end, score = line.split()            
                strand = "+"
                start, end = int(start), int(end)
                intron_start, intron_end = int(intron_start), int(intron_end)

                if start > end:
                    strand = "-"
                    start, end = end, start
                    intron_start, intron_end = intron_end, intron_start
            
                description = "tRNA-%s" %aa
                tno += 1
                id = "%s%04d" %(prefix, tno)
                size = end - start + 1 - (intron_end - intron_start + 1)
                        
                if intron_start != 0 :
                    pos = [start, intron_start, intron_end, end]
                else:
                    pos = [start, end]

            records += gff_record(id, scfid, source, feature, pos, score, strand, description, size)
    return records


parser = argparse.ArgumentParser(description='Convert tRNA/rRNA output to gff file')
parser.add_argument('result_file')
parser.add_argument('gff_outfile')
parser.add_argument('--source', dest='source', default = "rnammer", type = str)
parser.add_argument('--prefix', dest='prefix', default = "GL50803_r2", type = str)
args = parser.parse_args()

result_file = args.result_file
gff_outfile = args.gff_outfile
source = args.source
prefix = args.prefix

# Global variables
rno = 0
tno = 0

# process_RNA_file(tRNA_file, "tRNAscan-SE", "tRNA")
# process_RNA_file(tRNA_file, "rnammer", "rRNA")
with open(gff_outfile, 'w') as outh:
    print >>outh, process_RNA_file(result_file, source, prefix)
