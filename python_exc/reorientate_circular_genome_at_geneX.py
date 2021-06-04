#!/usr/bin/env python3
'''
    Re-orientate circular genome to start at 100bp upstream of specified gene, like dnaA, repB 
    or at the specified sequence 
    
'''
import argparse
from Bio import SeqIO


def get_gene(f):
    if f:
        gene = f.qualifiers.get('gene')[0] if "gene" in f.qualifiers else ''
        if gene.find('_') >= 0:
            gene = '_'.join(gene.split('_')[0:-1])
        if not gene:
            gene = f.qualifiers.get('product')[0] if "product" in f.qualifiers else ''
        return gene
    else:
        return ''

def get_pos(f):
    if f:
        start = f.location.start + 1 
        end = f.location.end
        strand = f.location.strand
        return start, end, strand
    else:
        return 0, 0, 0

parser = argparse.ArgumentParser(description='Re-orient fasta according to provided start gene or sequence ')
parser.add_argument('infile', help='Fasta file to re-orient')
parser.add_argument('--start_at', dest='start_at', type=str, help='genename or sequence to re-orient')
parser.add_argument('--gbk_file', dest='gbk_file', help='annotation file in genbank format')
parser.add_argument('--extra_bp', dest='extra_bp', default = 100, type=int, help='Extra bp to include for the start')


args = parser.parse_args()

fasta_file = args.infile
gbk_file = args.gbk_file
start_at = args.start_at
extra_bp = args.extra_bp
outfile = fasta_file + '.reoriented'


seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

def re_orient_to_start(start_at):
    orient_by_gene = True
    allowed = set('atcgn')
    if len(start_at) > 10 and set(start_at.lower()) <= allowed:
        orient_by_gene = False
    
    re_oriented = {}
    # print(start_at, orient_by_gene)
    if orient_by_gene:
        # Find out the pos of the given gene
        with open(gbk_file, 'r') as gbk_inh:
            for record in SeqIO.parse(gbk_inh, 'genbank'):
                seqid = record.id
                # print(seqid)
                for f in record.features:
                    if f.type != 'CDS':
                        continue
                    gene = get_gene(f)
                    if gene.startswith(start_at):
                        start, end, strand = get_pos(f)
                        seq = seq_dict[seqid]
                        if strand == -1:
                            # Make reverse complement
                            new_seq = seq[: end + extra_bp].reverse_complement() + seq[end + extra_bp :].reverse_complement()
                        else:
                            new_seq = seq[start - 1 - extra_bp :] + seq[: start - 1 - extra_bp]
                        print(seqid, start, strand)
                        re_oriented[seqid] = new_seq
    else:
        # Find out the pos of the given seq
        print('Im here')
        with open(fasta_file, 'r') as inh:
            for record in SeqIO.parse(inh, 'fasta'):
                seqid = record.id
                seq = record.seq.upper()
                start = seq.find(start_at.upper())
                if start < 0:
                    # Look for match on reverse complement strand
                    seq = seq.reverse_complement()
                    start = seq.find(start_at.upper())
                if start > 0:
                    # matches to the given seq
                    new_seq = seq[start :] + seq[: start]
                    print(seqid, start)
                    re_oriented[seqid] = new_seq
    return(re_oriented)



re_oriented = re_orient_to_start(start_at)

with open(fasta_file, 'r') as inh, open(outfile, 'w') as outh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id.split("|")[0]
        outh.write(">%s\n" %seqid)
        if seqid in re_oriented.keys():
            new_seq = re_oriented[seqid]
            outh.write(str(new_seq))
            print(len(new_seq))
        else:
            outh.write(str(record.seq))
        outh.write("\n")
        
