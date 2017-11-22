#!/usr/bin/env python

''' Convert Blast to gff3 output
    Made for rRNA.gff
    # Remove replicates
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Convert blast tab to gff file')
parser.add_argument('blast_tab')
parser.add_argument('scf_file')
parser.add_argument('gff_outfile')
parser.add_argument('seq_outfile')
parser.add_argument('--feature', dest='feature', default = "rRNA", type = str)
parser.add_argument('--prefix', dest='prefix', default = "GL50803_r2", type = str)
args = parser.parse_args()

blast_file = args.blast_tab
scf_file = args.scf_file
gff_outfile = args.gff_outfile
seq_outfile = args.seq_outfile
feature = args.feature
prefix = args.prefix

ref_seqs = defaultdict(Seq)
with open(scf_file, 'r') as fh:
    for record in SeqIO.parse(fh, 'fasta'):
        ref_seqs[record.id] = record.seq


def gff3_record(scfid, start, end, strand, score, evalue, description, id, partial, source = "BLAST", type="rRNA"):
    ''' Make a RNA gff3 record'''
    
    size = end - start + 1
    gene_attributes = "ID=%s;size=%d;description=%s" %(id, size, description)
    RNA_id = "rna_%s-1" %id
    RNA_attributes = "ID=%s;Parent=%s;size=%d;description=%s" %(RNA_id, id, size, description)
    
       
    exon_id = "exon_%s-1" %id
    exon_attributes = "ID=%s;Parent=%s;description=exon" %(exon_id, RNA_id)
    record =  "\t".join(map(str, [scfid, source, "gene", start, end, score, strand, ".", gene_attributes])) + "\n"
    record += "\t".join(map(str, [scfid, source, type, start, end, ".", strand, ".", RNA_attributes])) + "\n"
    record += "\t".join(map(str, [scfid, source, "exon", start, end, ".", strand, ".", exon_attributes]))
    return record
    


num = 0
na_records = []
done = defaultdict(list)
with open(blast_file, 'r') as fh, open(gff_outfile, 'w') as outh, open(seq_outfile, 'w') as seq_fh:
    for line in fh:
        partial = None
        arr = line.strip().split("\t")
        qseqid, sseqid = arr[0:2]
        qstart, qend, sstart, send = map(int, arr[6:10])
        evalue, score = map(float, arr[10:12])
        stitle = arr[-1]
        slen = int(arr[13])
        seq_len = len(ref_seqs[qseqid])
        if sstart > send:
            strand = "-"
            sstart, send = send, sstart
        else:
            strand = "+"         
        
        if strand == "+":
            if sstart != 1:
                qstart -= sstart -   1
                if qstart < 1:
                    qstart = 1
                    partial = 5
            if send != slen:
                qend += slen - send
                if qend > seq_len:
                    qend = seq_len
                    partial = 3
            na_seq = ref_seqs[qseqid][qstart-1:qend]
        else:
            if sstart != 1:
                qend += (sstart - 1)
                if qend > seq_len:
                    qend = seq_len
                    partial = 5
            if send != slen:
                qstart -= slen - send
                if qstart < 1:
                    qstart = 1
                    partial = 3
        
            na_seq = ref_seqs[qseqid][qstart-1:qend].reverse_complement()
        
        # # double check 5S end seq, and adjust accordlingly
        # if stitle.startswith("5S"):
        #     if not na_seq.endswith("ATTTT"):
        #         index = na_seq.find("ATTTT")
        #         if index < 0 or index > len(na_seq) - 5:
        #             print line
        #         if strand == "+":
        #             qend = qstart + index + 4
        #             na_seq = ref_seqs[qseqid][qstart-1:qend]
        #         else:
        #             qstart = qend - (index + 4)
        #             na_seq = ref_seqs[qseqid][qstart-1:qend].reverse_complement()
        

        if (qseqid, strand) not in done:
            done[(qseqid, strand)].append([qstart, qend])
        elif [qstart, qend] in done[(qseqid, strand)]:
            continue
        else:
            done[(qseqid, strand)].append([qstart, qend])
                    
        num += 1
        id = "%s%04d" %(prefix, num)
        print >>outh, gff3_record(qseqid, qstart, qend, strand, score, evalue, stitle, id, partial, "Blast", feature)
    
        na_record = SeqRecord(Seq(str(na_seq), IUPAC.ambiguous_dna), id = id, description = stitle)
        na_records.append(na_record)
        
        
    SeqIO.write(na_records, seq_fh, 'fasta')



    
    
    

        
            