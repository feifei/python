#!/usr/bin/env python

''' Convert Blast to gff3 output
    # Remove replicates
'''
import os 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Convert blast tab to gff file')
parser.add_argument('blast_tab')
parser.add_argument('--feature', dest='feature', default = "gene", type = str)
parser.add_argument('--prefix', dest='prefix', default = "X", type = str)
parser.add_argument('--evalue', dest='evalue', default = 1e-10, type = float)
parser.add_argument('--color', dest='color', default = "2", type = str)
args = parser.parse_args()

blast_tab = args.blast_tab
feature = args.feature
prefix = args.prefix
evalue_cutoff = args.evalue
color = args.color

basename, extension = os.path.splitext(blast_tab)

gff_file = basename + ".gff"

num = 0
with open(gff_file, 'w') as gff_h, open(blast_tab, 'r') as tab_h:
    for line in tab_h:
        arr = line.strip().split("\t")
        qid, sid = arr[0:2]
        identity = arr[2]
        aln_len = int(arr[3])
        qstart, qend, sstart, send = map(int, arr[6:10])
        evalue = float(arr[10])
        bitscore = arr[11]
        
        if evalue > evalue_cutoff:
            continue
            
        # if reverse:
        #     seqid, target = target, seqid
        #     qstart, sstart = sstart, qstart
        #     qend, send = send, qend
        
        num += 1
        ID = "%s-%d" %(prefix, num)
        attribs = "ID=%s;Target=%s;evalue=%.1e;aln_len=%d;identity=%s;qry_start=%d;qry_end=%d;" %(ID, qid, evalue, aln_len, identity, qstart, qend)
        
        if color:
            attribs += "color=%s" %(color)


        if sstart < send:
            strand = "+"
        else:
            strand = "-"
            sstart, send = send, sstart
            
        
        print >>gff_h, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t.\t%s" %(sid, "BLAST", feature, sstart, send, bitscore, strand, attribs)