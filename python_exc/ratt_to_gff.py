#!/usr/bin/env python
'''
    Convert RATT output to a single GFF3 file, adressing
    some problems regarding RATT transferred annotation, including
    - Genes doesn't start with a start codon M
    - Genes contain splice sites
    - Genes short (Maybe comparison with the origin length)
    - Genes overlap each other
    - Introns were handled manually
'''
import os, re
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO

from overlap import *

parser = argparse.ArgumentParser(description='Convert ratt run to one gff file')
parser.add_argument('ratt_run_dir')
parser.add_argument('scf_file')
parser.add_argument('gff_file')
parser.add_argument('--gcode', dest='genetic_code', default = 1, type = int)
args = parser.parse_args()

run_dir = args.ratt_run_dir
gff_file = args.gff_file
scf_file = args.scf_file
genetic_code = args.genetic_code

seq_dict = defaultdict(Seq)
with open(scf_file, 'r') as scf_fh:
    for record in SeqIO.parse(scf_fh, 'fasta'):
        seq_dict[record.id] = record.seq


def parse_position(p):
    if re.match("complement", p):
        strand = "-"
    else:
        strand = "+"
    start, end = re.search("(\d+)\.\.(\d+)", p).groups()
    return (int(start), int(end), strand)


def get_aa_seq(seqid, strand, start, end, genetic_code=genetic_code):
    gene_seq = seq_dict[seqid][start - 1 : end]
    if strand == "-":
        gene_seq = gene_seq.reverse_complement()
    return gene_seq.translate(genetic_code)
    

def make_gff(seqid, start, end, strand, geneid, description):
    attributes = "ID=%s;description=%s;label=%s;color=3" %(geneid, description, geneid.split("_")[1] + "_" + description)
    return "\t".join(map(str, [seqid, "RATT", "gene", start, end, ".", strand, ".", attributes]))


strand = ""
d = defaultdict(list)
with open(gff_file, 'w') as gff_fh:
    for dirpath, _, filenames in os.walk(run_dir):
        for filename in filenames:
            if filename.endswith("final.embl"):
                scfid = filename.split(".")[1]
                infile = os.path.join(dirpath, filename) 
                with open(infile, 'r') as inh:
                    # EMBL format imperfect, so parse manually instead of SeqIO
                    for line in inh:
                        line = line.rstrip()
                        if not line.startswith("FT"):
                            continue
                        if line.startswith("FT   CDS"):
                            pos = line.split()[2]
                            if re.search("join", pos):
                                # fragmented
                                inh.next()
                                inh.next()
                                inh.next()
                                continue

                            start, end, strand = parse_position(pos)
                            aa_seq = get_aa_seq(scfid, strand, start, end)
                    
                    
                        if re.search("\/gene=", line):
                            while not line.endswith("\""):
                                line = line + " " + inh.next().split("FT")[1].strip()
                            gene = line.split("\"")[1]
                            description, geneid = re.match("(.*?) \[(.*?)\]", gene).groups()
                                                
                            # Ignore tRNA, rRNA, will be done separately
                            if geneid.split("_")[1].startswith("r") or geneid.split("_")[1].startswith("t"):
                                continue
                        
                            if re.search("deprecated", description):
                                continue

                            if [start, end] not in d[(scfid, strand)]:
                                d[(scfid, strand)].append([start, end])
                            else:
                                # Completely overlap, skip
                                continue

                            if not aa_seq.startswith("M"):
                                print geneid, scfid, "NOT starting with M", strand, start, end, description
                            elif not aa_seq.endswith("*"):
                                print geneid, scfid, "NOT ending with *", strand, start, end, description
                            elif aa_seq.strip("*").find("*") != -1:
                                print geneid, scfid, "contains *", strand, start, end, description
                            gff_str = make_gff(scfid, start, end, strand, geneid, description)
                            print >>gff_fh, gff_str
                          
                    
                            
                                