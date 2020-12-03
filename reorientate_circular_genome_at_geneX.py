#!/usr/bin/env python
'''
    Re-orientate circular genome to start at dnaA + 100bp, like CP048439 genome
    
'''
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', help="Fasta file to re-orient")
parser.add_argument('blast_file', help="blast tab result")
args = parser.parse_args()

fasta_file = args.infile
outfile = fasta_file + '.reoriented'
blast_file = args.blast_file

blast_inh = open(blast_file, 'r')
arr = blast_inh.readline().split("\t")
qseqid, sseqid = arr[0:2] 
qstart, qend, sstart, send = map(int, arr[6:10])
qstart, qend
blast_inh.close()        


def re_oriente_to_dnaA(ref_seq, start, end, sstart, ssend, extra = 100):
    if sstart < ssend:
        reverse = False
    else:
        reverse = True
    
    if reverse:
        # ref_seq should be reverse complemeneted
        return ref_seq[:end+100].reverse_complement() + ref_seq[end+100:].reverse_complement()
        
    else:
        return ref_seq[start-1-100:] + ref_seq[:start-1-100]


with open(fasta_file, 'r') as inh, open(outfile, 'w') as outh:
    for record in SeqIO.parse(inh, 'fasta'):
        seqid = record.id.split("|")[0]
        outh.write(">%s\n" %seqid)
        if seqid == qseqid:
            new_seq = re_oriente_to_dnaA(record, qstart, qend, sstart, send)
            outh.write(str(new_seq.seq))
            print(len(new_seq.seq))
        else:
            outh.write(str(record.seq))
        outh.write("\n")
        
