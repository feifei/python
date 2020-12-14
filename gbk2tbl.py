#!/usr/bin/env python3
'''
    Convert gbk file from prokka annotation to .tbl in preparation for submission to NCBI
'''

import re, os
import argparse
from collections import defaultdict
from Bio import SeqIO


def prepare_repeat_feature(start, end, strand, t, note="", rpt_family="", rpt_type="", rpt_unit_seq="", pseudos=[]):
    ''' Prepare repeat feature '''
    
    def feature_head(start, end, tag):
        return "%d\t%d\t%s\n" %(start, end, tag)

    def feature_line(tag, tag_content):
        return "\t\t\t%s\t%s\n" %(tag, tag_content)


    return feature_head(start, end, t) + feature_line("note", note) + \
    feature_line('rpt_family', rpt_family) + feature_line('rpt_type', rpt_type) + \
    feature_line('rpt_unit_seq', rpt_unit_seq)
    
    
    
def prepare_feature(geneid, start, end, strand, product, t, gene = "", dbname = "FXUU"):
    ''' Prepare gene CDS feature for each CDS '''
    
    def feature_head(start, end, tag):
        return "%d\t%d\t%s\n" %(start, end, tag)

    def feature_line(tag, tag_content):
        return "\t\t\t%s\t%s\n" %(tag, tag_content)
    
    def feature_without_intron(t, start, end, product, protein_id, transcript_id):
        return feature_head(start, end, t) + feature_line("product", product) + \
               feature_line("protein_id", protein_id) + feature_line("transcript_id", transcript_id) 
    
    def add_pseudo():
        return "\t\t\tpseudo\n" # + feature_line("note", "nonfunctional due to frameshift")
    

    start, end = (end, start) if strand == -1 else (start, end)
    protein_id = "gnl|%s|%s" %(dbname, geneid)
    if gene:
        gene_feature = feature_head(start, end, "gene") + feature_line("gene", gene) + feature_line("locus_tag", geneid)
    else:
        gene_feature = feature_head(start, end, "gene") + feature_line("locus_tag", geneid)

    transcript_id = "gnl|%s|%s.%s" %(dbname, t.lower(), geneid)
    
    if t == "CDS": 
        cds_feature = feature_without_intron("CDS", start, end, product, protein_id, transcript_id)
    else:
        cds_feature = feature_without_intron(t, start, end, product, protein_id, transcript_id)
    
    if geneid in pseudos: # H8529_04859
        #pseudo
        return gene_feature + feature_line('gene_desc', product) + add_pseudo()
    else:    
        #Bacteria has no mrna
        return gene_feature + cds_feature 



parser = argparse.ArgumentParser(description='Convert .gbk to .tbl for submission to NCBI')
parser.add_argument('infile', help='gbk file to convert')
parser.add_argument('--pseudos', dest='pseudos', default='', type=str, help='pseudo geneids, seprated by ,')
parser.add_argument('--seqids', dest='seqids', default='', help='sequence ids to specify, seprated by ,')
parser.add_argument('--locus_tag_prefix', dest='locus_tag_prefix', default='', help='Locus tag prefix')


args = parser.parse_args()

gbk_file = args.infile
pre, ext = os.path.splitext(gbk_file)
tbl_file = pre + '.tbl'
pseudos = args.pseudos.split(',')
seqids = args.seqids.split(',')
locus_tag_prefix = args.locus_tag_prefix

#ids = ['EC93', 'pCP127']
i = 0
rna_start, rna_end, cds_start, cds_end = 0, 0, 0, 0
prev_feature_to_output = ""
with open(gbk_file, 'r') as inh, open(tbl_file, 'w') as outh:
    for record in SeqIO.parse(inh, "genbank"):
        seqid = record.id.split('_genome')[0]
        if seqids != ['']:
            print(seqids)
            seqid = seqids[i]
            i += 1
        for f in record.features:
            if f.type == 'source':
                outh.write(prev_feature_to_output)
                prev_feature_to_output = ""
                outh.write(">Features %s\n" %seqid)
                continue
            skip_tag=False
            prev_skip_tag = False
            start = f.location.start + 1 
            end = f.location.end
            strand = f.location.strand
            # print(start, end, strand)
            t = f.type 
            
            if t == 'CDS':
                cds_start = start
                cds_end = end
                if cds_start >= rna_start and cds_end <= rna_end:
                    # CDS within misc_RNA, should ski misc_RNA
                    print('CDS within misc_RNA', f.qualifiers.get('locus_tag'), f.qualifiers.get('gene'))
                    prev_skip_tag = True
            elif t == 'misc_RNA':
                rna_start = start
                rna_end = end
                if rna_start >= cds_start and rna_end <= cds_end:
                    # misc_RNA within CDS skip
                    print('Skipped, misc_RNA within CDS', f.qualifiers.get('locus_tag'), f.qualifiers.get('gene'))
                    skip_tag = True
                    
            if t == 'repeat_region':
                note = f.qualifiers.get('note').pop()
                rpt_family = f.qualifiers.get('rpt_family').pop()
                rpt_type = f.qualifiers.get('rpt_type').pop()
                rpt_unit_seq = f.qualifiers.get('rpt_unit_seq').pop()
                outh.write(prepare_repeat_feature(start, end, strand, t, note, rpt_family, rpt_type, rpt_unit_seq))
            else:
                if skip_tag:
                    # skip misc_RNA within CDS 
                    continue
                locus_tag = f.qualifiers.get('locus_tag').pop() if 'locus_tag' in f.qualifiers else ''
                if locus_tag_prefix:
                    locus_tag = locus_tag_prefix + '_' + locus_tag.split('_')[1]
                product = f.qualifiers.get('product').pop() if 'product' in f.qualifiers else ''
                gene = f.qualifiers.get('gene').pop() if "gene" in f.qualifiers else ''
                                    
                feature_to_output = prepare_feature(locus_tag, start, end, strand, product, t, gene)
                if prev_skip_tag or prev_feature_to_output == "":
                    # Skip CDS within misc_RNA
                    prev_feature_to_output = feature_to_output
                    continue
                
                outh.write(prev_feature_to_output)
                prev_feature_to_output = feature_to_output

    outh.write(prev_feature_to_output)            