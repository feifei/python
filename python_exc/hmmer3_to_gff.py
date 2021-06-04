#!/usr/bin/env python

''' Convert hmmer3 output to gff format'''

from collections import defaultdict
import argparse
import os, re

from HMMER import *
from handy import *

parser = argparse.ArgumentParser()
parser.add_argument('genes_gff', help = "genes to annotate in gff file")
parser.add_argument('pfam_out', help = "Pfam results to annotate the genes domains")
args = parser.parse_args()

gene_file = args.genes_gff
pfam_file = args.pfam_out

basename, extension = os.path.splitext(pfam_file)
pfam_gff_file = basename + ".gff"


genes_info = defaultdict(list)
with open(gene_file, 'r') as gene_fh:
    for line in gene_fh:
        line = line.strip()
        if len(line) == 0:
            continue
        arr = line.split("\t")
        scfid = arr[0]
        start, end = map(int, arr[3:5])
        strand = arr[6]
        geneid = re.match("ID=(.*?);", arr[-1]).groups()[0]
        genes_info[geneid] = [scfid, start, end, strand]
    


def gff3_record(scfid, geneid, begin, end, score, evalue, strand, target, description, targe_begin, target_end, color="", source="HMMER3", typ="match"):
    ''' Make a gff3 record '''    
    
    gene_attributes = "ID=%s;Target=%s;Description=%s;score=%.2f;evalue=%s;color=%s;" %(geneid, target, description.strip(), score, str(evalue), color)    
    record =  "\t".join(map(str, [scfid, source, typ, begin, end, ".", strand, ".", gene_attributes]))
    return record


good_total = 0
weak_total = 0

# EVAL_BIG_CUTOFF = 0.01
# EVAL_SAMLL_CUTOFF = 1e-05
SCORE_SMALL_CUTOFF = 15
SCORE_BIG_CUTOFF = 20
WEAK_COLOR = "255 211 155"
GOOD_COLOR = "244 164 96"

with open(pfam_gff_file, 'w') as gff_h, open(pfam_file, 'r') as hmmer_h:
    for hmmresult in parseMultiHMMER3(hmmer_h):
        geneid = hmmresult.seqName
        scfid, scf_begin, scf_end, strand = genes_info[geneid]
        
        last_name = ""
        for hmmunit in hmmresult.units:
            na_begin, na_end = aa_to_na_pos(hmmunit.seqFrom, hmmunit.seqTo)
            begin, end = gene_to_scf_pos(strand, scf_begin, scf_end, na_begin, na_end)            
            
            score = hmmunit.bits
            # dom_evalue = hmmunit.domEvalue
            # domain = hmmunit.domain
            evalue = hmmunit.evalue
            # accuracy = hmmunit.aliAcc
            name = hmmunit.name
            description = hmmresult.seqs[name].desc
            target_begin = hmmunit.hmmFrom
            target_end = hmmunit.hmmTo
            if name != last_name and last_name != "":
                print >>gff_h, ""
            last_name = name
        
            # Setting criteria, individual domain score >= 15
            if score >= SCORE_BIG_CUTOFF:
                print >>gff_h, gff3_record(scfid, geneid, begin, end, score, evalue, strand, name, description, target_begin, target_end, GOOD_COLOR, "HMMER3_good")
                good_total += 1
            elif score >= SCORE_SMALL_CUTOFF:
                print >>gff_h, gff3_record(scfid, geneid, begin, end, score, evalue, strand, name, description, target_begin, target_end, WEAK_COLOR, "HMMER3_weak")
                weak_total += 1
                    
                
print "There are %d good hits, and %d weak hits in total" % (good_total, weak_total)

